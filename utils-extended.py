SAME_G = False
func = lambda x: np.sqrt(x)
func2 = lambda x: np.sqrt(x)

def F(sol, obs_util, m):
    v = 0

    for j in range(m):
        tmp = np.sum(sorted([obs_util[i][j] for i in sol], reverse=True))
        if j != 1:
            v += func(tmp) * weight_F[j]
        else:
            if not SAME_G:
                v += func2(tmp) * weight_F[j]
            else:
                v += func(tmp) * weight_F[j]
    return v

def marg_F(sol, item, obs_util, m):
    sol2 = copy.deepcopy(sol)
    sol2.append(item)
    return F(sol2, obs_util, m) - F(sol, obs_util, m)


def marg_linear_mult(sol, items, obs_util, m, pr=False):
    v = np.array([0 for j in range(m)])

    for j in range(m):
        tmp = [obs_util[i][j] for i in sol]
        v[j] += np.sum(tmp)

    cur_util = 0
    for j in range(m):
        cur_util += v[j] * weight_F[j]

    utils = []
    for i in items:
        if i in sol:
            utils.append(0)
        else:
            assert(obs_util[i][j] >= 0)
            new_util = 0
            for j in range(m):
                new_util += (v[j] + obs_util[i][j]) * weight_F[j]

            assert(new_util >= cur_util)
            utils.append(new_util - cur_util)
    return utils


def marg_F_mult(sol, items, obs_util, m, pr=False):
    v = np.array([0 for j in range(m)])

    for j in range(m):
        tmp = [obs_util[i][j] for i in sol]
        v[j] += np.sum(tmp)

    cur_util = 0
    for j in range(m):
        if j != 1:
            cur_util += func(v[j]) * weight_F[j]
        else:
            if not SAME_G:
                cur_util += func2(v[j]) * weight_F[j]
            else:
                cur_util += func(v[j]) * weight_F[j]


    utils = []
    for i in items:
        if i in sol:
            utils.append(0)
        else:
            assert(obs_util[i][j] >= 0)
            new_util = 0
            for j in range(m):
                if j != 1:
                    new_util += func(v[j] + obs_util[i][j]) * weight_F[j]
                else:
                    if not SAME_G:
                        new_util += func2(v[j] + obs_util[i][j]) * weight_F[j]
                    else:
                        new_util += func(v[j] + obs_util[i][j]) * weight_F[j]
            assert(new_util >= cur_util)
            utils.append(new_util - cur_util)
    return utils

def get_top_n(items, utils, k):
    pairs = list(zip(utils, items))
    sol = sorted(pairs, reverse=True)[:k]
    return [s[1] for s in sol]

def get_top_n_set(items, utils, k):
    return set(get_top_n(items, utils, k))

DEBUG = False

def print_details(item, lat_util, obs_util):
    print(f'artist#{item} {np.round(lat_util[item][0],2)},\t{lat_util[item][1]},\t{lat_util[item][2]})')
    print(f'artist#{item} {np.round(obs_util[item][0],2)},\t{obs_util[item][1]},\t{obs_util[item][2]})')
