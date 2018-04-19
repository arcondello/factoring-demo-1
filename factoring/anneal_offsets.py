"""
Utilities for finding anneal offsets for chains.
"""

from __future__ import division

from six import iteritems
import numpy as np
import scipy.optimize as opt


def build_chain_offsets(anneal_offset_ranges, embedding, c=2):
    """
    Builds a set of anneal offsets for an embedding, based on the strategy
    defined in https://www.dwavesys.com/sites/default/files/14-1002A_B_tr_Boosting_integer_factorization_via_quantum_annealing_offsets.pdf

    For this function, the annealing offset of a variable is a function
    of its chain length. The idea is that longer chains freeze out
    earlier, so they should be given an anneal delay. The delay on a chain
    of length k is proportional to c^((k-1)/k)-1 for some constant c >= 1.

    Args:
        anneal_offset_ranges (list of lists):
            available offset ranges on each qubit, from
            solver.properties.anneal_offset_ranges. offset_ranges[i] has the
            form [min_offset, max_offset] where min_offset <= 0, max_offset
            >= 0.

        embedding (dict):
            an embedding, i.e. a dictionary mapping variables to chains of
            qubits (for example as returned by minorminer.find_embedding).
            Unused qubits are assigned an anneal offset of 0.

        c (float, default = 2):
            The base in the delay function delay(k) = c^((k-1)/k)-1. c
            should be >= 1. Experimentally, c = 2 or 1/log 2 are reasonable
            choices.

    Returns:
        offset_strategy (list):
            a list of the anneal offsets for each qubit. The output of the
            function can be passed into a sampler using the
            optional 'anneal_offsets' argument.

    Given the working offset strategy that is output here, you will want
    to sweep over the scale of the offset effect:

        working_strategy = working_offset_strategy(...)
        for i in range(11):
            sample_args['anneal_offsets'] = (i/10.)*working_strategy
            answer[i] = sample_ising(...,sample_args)
        end
        pick best answer
    """

    ideal_offsets = chain_offsets_nonlinear(embedding, c)
    return working_offset_strategy(anneal_offset_ranges, ideal_offsets)


def embed_offsets(embedding, var_offsets):
    """
    Given an offset function for variables, write it as an
    offset function for qubits.

    Args:
        embedding (dict):
            an embedding, i.e. a dictionary mapping variables to chains of
            qubits (for example as returned by minorminer.find_embedding).

        var_offsets (dict):
            a dictionary mapping each variable to an anneal offset.

    Returns:
        qubit_offsets (dict):
            a dictionary mapping qubits to anneal offsets. Only qubits used
            in the embedding are included.
    """
    qubit_offsets = dict()
    for v, offset in iteritems(var_offsets):
        if v in embedding:
            for s in embedding[v]:
                qubit_offsets[s] = offset

    return qubit_offsets


def chain_offsets_nonlinear(embedding, c=2):
    """
    Determine a set of anneal offsets for variables in an embedding,
    based on the strategy defined in based on the strategy
    defined in https://www.dwavesys.com/sites/default/files/14-1002A_B_tr_Boosting_integer_factorization_via_quantum_annealing_offsets.pdf

    For this function, the annealing offset of a variable is a function
    of its chain length. The idea is that longer chains freeze out
    earlier, so they should be given an anneal delay. The delay on a chain
    of length k is proportional to c^((k-1)/k)-1 for some constant c >= 1.

    The anneal offset is normalized to the range [-1,0]. Positive output is
    advance, negative output is delay. A chain of length 1 should have
    offset 0, and a chain of maximum length should have offset -1 (after
    normalization).

    Args:
        embedding (dict):
            an embedding, i.e. a dictionary mapping variables to chains of
            qubits (for example as returned by minorminer.find_embedding).

        c (float, default = 2):
            The base in the delay function delay(k) = c^((k-1)/k)-1. c
            should be >= 1. Experimentally, c = 2 or 1/log 2 are reasonable
            choices.

    Returns:
        ideal_offsets (dict):
            a dictionary mapping qubits to anneal offsets. The output of
            this function can be passed into working_offset_strategy().
    """

    def chain_delay(x, c):
        # delay for chains of length x (c is base)
        return c ** ((x - 1) / x) - 1

    delay = dict()
    for v, chain in iteritems(embedding):
        x = len(chain)
        if x > 0:
            delay[v] = chain_delay(x, c)

    # normalize to [0,1], and define advance as -delay:
    min_delay = 0
    max_delay = max(delay.values())
    advance = {k: -(v-min_delay)/(max_delay-min_delay)
               for k, v in iteritems(delay)}

    # translate variable offsets into qubit offsets:
    qubit_offsets = embed_offsets(embedding, advance)

    return qubit_offsets


def working_offset_strategy(anneal_offset_ranges, ideal_offsets,
                            strategy_type='full'):
    """
    From an ideal offset strategy in the range [-1,1] (for example the
    output of chain_offsets_nonlinear), scale and shift to produce a working
    offset strategy using the available offset ranges in the solver such that
    the scale of the offsets is as large as possible.

    Args:
        anneal_offset_ranges (list of lists):
            available offset ranges on each qubit, from
            solver.properties.anneal_offset_ranges. offset_ranges[i] has the
            form [min_offset, max_offset] where min_offset <= 0, max_offset
            >= 0.

        ideal_offsets (dict):
            a dictionary mapping qubits to anneal offsets (for example
            from chain_offsets_nonlinear).

        strategy_type (string):
            either:
                "full" = use all available offset range. Requires scipy.
                "zeroed" = use all available offset range, subject to 0
                           being centered at 0.
                In general, if you shift the offset on every single qubit by a
                constant, then it won't have any effect. So "full" should be
                strictly better than "zeroed". However "zeroed" does not
                require scipy.

    Returns:
        offset_strategy (list):
            a list of the anneal offsets for each qubit. The output of the
            function can be passed into a sampler using the
            optional 'anneal_offsets' argument.


    Given the working offset strategy that is output here, you will want
    to sweep over the scale of the offset effect:

        working_strategy = working_offset_strategy(...)
        for i in range(11):
            sample_args['anneal_offsets'] = (i/10.)*working_strategy
            answer[i] = sample_ising(...,sample_args)
        end
        pick best answer
    """

    # convert offset_ranges to numpy (allows access by list of indices)
    offset_ranges = np.array(anneal_offset_ranges)

    if not ideal_offsets:
        scale = 0
        shift = 0
    elif strategy_type == 'full':
        # strategy "full": use as much range as possible.
        # maximize "scale" subject to scale*offset + shift in range.
        # This is an LP in the variables x = (scale,shift).

        used_qubits = ideal_offsets.keys()

        # Ax <= b: each row corresponds to a qubit and looks like
        # [ideal_offset, 1]*[scale, shift]' <= maximum_offset.
        A1 = np.ones(shape=(len(used_qubits), 2))
        A1[:, 0] = np.array(ideal_offsets.values())
        b1 = offset_ranges[used_qubits, 1]

        # similarly for [ideal_offset, 1]*[scale, shift]' >= minimum_offset.
        A2 = -A1
        b2 = -offset_ranges[used_qubits, 0]

        # concatenate inequalities
        A_ub = np.concatenate((A1, A2), axis=0)
        b_ub = np.concatenate((b1, b2))

        # objective function: maximize scale
        # (i.e. minimize [-1,0]*[scale, shift]'.)
        c = np.array([-1, 0])

        # bounds (0 <= scale <= Inf, -Inf <= shift <= Inf):
        bounds = [(0, None), (None, None)]

        # solve (uses scipy linear programming):
        result = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds)

        # extract result.
        if result.status != 0:
            raise ValueError(result.message)
        x = result.x
        scale = x[0]
        shift = x[1]

    elif strategy_type == 'zeroed':
        # strategy "zeroed": use as much range as possible, subject to 0
        # offset being centered at 0. That is, maximize "scale" subject
        # to "shift" = 0.

        # initialize the largest allowable positive and negative scales:
        cjj_problem_max = np.Inf
        cjj_problem_min = np.Inf

        for k, v in iteritems(ideal_offsets):
            if v < 0:
                # update the negative scale:
                cjj_problem_min = min(cjj_problem_min, offset_ranges[k, 0]/v)
            elif v > 0:
                # update the positive scale:
                cjj_problem_max = min(cjj_problem_max, offset_ranges[k, 1]/v)

        # allowable scale the smallest of the negative and positive scales:
        scale = min(cjj_problem_max, cjj_problem_min)
        shift = 0

    else:
        raise ValueError('unrecognized strategy type')

    # define working offsets
    offset_strategy = {k: scale*v + shift for k, v in iteritems(ideal_offsets)}

    # fix rounding/truncating errors:
    offset_strategy = {k: min(max(v, offset_ranges[k, 0]), offset_ranges[k, 1])
                       for k, v in iteritems(offset_strategy)}

    # convert to API format (a list of offsets, one for each qubit):
    offset_list = [0]*len(anneal_offset_ranges)
    for k, v in iteritems(offset_strategy):
        offset_list[k] = v

    return offset_list
