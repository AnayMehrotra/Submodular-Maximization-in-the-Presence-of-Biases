def baseline_greedy_uncons(obs_util, grps, k, m, pr=False):
    n = np.sum([len(g) for g in grps])
    all = [i for i in range(n)]

    sol = [] #()

    for i in range(k):
        utils = marg_F_mult(sol, all, obs_util, m)
        tmp_sol = get_top_n(all, utils, 1)

        assert(len(tmp_sol) == 1)
        ind = list(tmp_sol)[0]
        all.remove(ind)

        sol.append(ind)

    assert(len(sol)==k)
    return sol

def baseline_greedy_rooney(obs_util, grps, k, m, pr=False):
    n = np.sum([len(g) for g in grps])
    all = [i for i in range(n)]

    sol = [] #()

    for i in range(k):
        utils = marg_F_mult(sol, all, obs_util, m)

        if i == k-1:
            disadv_selected = False
            for s in sol:
                if s in grps[1]:
                    disadv_selected = True
                    break

            if not disadv_selected:
                all = list(grps[1])

        tmp_sol = get_top_n(all, utils, 1)

        assert(len(tmp_sol) == 1)
        ind = list(tmp_sol)[0]
        all.remove(ind)
        sol.append(ind)

    assert(len(sol)==k)
    return sol

def baseline_greedy_equal(obs_util, grps, k, m):
    p = len(grps)
    n = np.sum([len(g) for g in grps])
    all = [i for i in range(n)]

    pr = False

    sol = [] #list()

    cnt = [0 for t in range(p)]

    while len(sol) < k:
        utils = marg_F_mult(sol, all, obs_util, m)
        tmp_sol = get_top_n(all, utils, 1)

        assert(len(tmp_sol) == 1)
        ind = list(tmp_sol)[0]

        inv_ind = -1
        for i in range(len(all)):
            if all[i] == ind:
                inv_ind = i
                break

        fg = 1
        for t in range(p):
            if ind in grps[t]:
                if cnt[t] > max(k * 0.5, k - len(grps[1-t])):
                    fg = 0

        if fg:
            sol.append(ind)
            for t in range(p):
                if ind in grps[t]:
                    cnt[t] += 1

        all.remove(ind)
    assert(len(sol)==k)
    return sol

def baseline_greedy_cons(obs_util, grps, k, m):
    p = len(grps)
    n = np.sum([len(g) for g in grps])
    all = [i for i in range(n)]

    pr = False

    sol = [] #list()

    cnt = [0 for t in range(p)]

    while len(sol) < k:
        utils = marg_F_mult(sol, all, obs_util, m)
        tmp_sol = get_top_n(all, utils, 1)

        assert(len(tmp_sol) == 1)
        ind = list(tmp_sol)[0]

        inv_ind = -1
        for i in range(len(all)):
            if all[i] == ind:
                inv_ind = i
                break

        fg = 1
        for t in range(p):
            if ind in grps[t]:
                if cnt[t] >= k * len(grps[t]) / n:
                    fg = 0

        if fg:
            sol.append(ind)
            for t in range(p):
                if ind in grps[t]:
                    cnt[t] += 1

        all.remove(ind)
    assert(len(sol)==k)
    return sol

def algo_3_disjoint_attr(obs_util, grps, k, m): # lat_util = None):
    p = len(grps)
    n = np.sum([len(g) for g in grps])
    grp_lists = [list(grps[t]) for t in range(p)]

    def count_sol_in_attr(sol, obs_util):
        k_pr = np.array([0.0 for j in range(m)])

        for j in range(m):
            for i in sol_tmp:
                k_pr[j] += (obs_util[i][j] > 0)

        k_pr *= float(float(k)/np.sum(k_pr))

        return k_pr

    #####################################

    obs_util_grp_a = copy.deepcopy(obs_util)
    for i in grps[1]: obs_util_grp_a[i] *= 0
    sol_tmp = baseline_greedy_uncons(obs_util_grp_a, grps, k * len(grps[0]) // n, m, pr=False)
    k_pr = count_sol_in_attr(sol_tmp, obs_util)

    ####################################

    sol = []

    for j in range(m):
        for t in range(p):
            k_int = int(k_pr[j] * len(grps[t]) / n)
            utils = [obs_util[i][j] for i in grps[t]]
            tmp_sol = get_top_n(grps[t], utils, k_int)
            for tmp in tmp_sol:
                sol.append(tmp)
    sol = set(sol)

    if len(sol) < k:
        unselec = []
        for i in range(n):
            if i not in sol:
                unselec.append(i)

        utils = marg_F_mult(sol, unselec, obs_util, m)
        tmp_sol = get_top_n(unselec, utils, k - len(sol))

        for i in tmp_sol:
            sol.add(i)
            if len(sol) == k:
                break
        sol = list(sol)

    assert(len(sol)==k)

    return sol
