def get_movie_dist_over_genre(u_id):
    
    rated_movies = get_rated_movie_indices(u_id)
    
    dist = np.zeros(len(genres))
    
    for m_ind in rated_movies:
        for j, g in enumerate(genres):
            if g not in genres_matching_tags: 
                continue
            dist[j] += movie_genres[m_ind][j]
    
    return dist / np.sum(dist)
    
def get_rated_movie_indices(u_id):
    rated_movies = []
    
    for rating in user_ratings_for_movies_with_scores[u_id]:
        if rating['movie_index'] in grps_set[0] or rating['movie_index'] in grps_set[1]:
            rated_movies.append(rating['movie_index'])
            assert(rating['movie_index'] in movie_ind_with_scores)
    
    return rated_movies

def get_score_rating_draw(sol, m_ids, m_ind_to_i, cnt=100):
    score = 0
    
    mx=-1
    mn=1e8
    
    user_draws = rng.choice(len(selected_users), cnt)
    
    cnt2 = 0
    
    
    for u_id in user_draws:
        if u_id not in user_ratings_for_movies_with_scores: continue 
        for rating in user_ratings_for_movies_with_scores[u_id]:
            m_ind = rating['movie_index'] 

            if m_ind not in m_ind_to_i: continue

            if m_ind_to_i[m_ind] in sol:
                score += rating['rating']
                cnt2 += 1
                mx = max(rating['rating'], mx)
                mn = min(rating['rating'], mn)
    
    score /= cnt2

    print(f'max: {mx}, min: {mn}')
    
    return score
    
    
def get_score_rating(sol, u_id, m_ids, m_ind_to_i):
    score = []
    
    mx=-1
    mn=1e8
    
    for rating in user_ratings_for_movies_with_scores[u_id]:
        m_ind = rating['movie_index'] 
        
        if m_ind not in m_ind_to_i: continue
        
        if m_ind_to_i[m_ind] in sol:
            score.append(rating['rating'])
            mx = max(rating['rating'], mx)
            mn = min(rating['rating'], mn)
    
    print(f'max: {mx}, min: {mn}')
    
    print(f'score-mean={np.mean(score)} and score-std={np.std(score)}')
    print(score)
    
    return np.sum(score)

def get_score_rating_overall_users(sol, m_inds, m_ind_to_i):
    score = []
    
    mx=-1
    mn=1e8

    for i in sol:
        m_ind = m_inds[i]
        score.append(movie_rating_sum[m_ind]/movie_rating_cnt[m_ind])
       
    return np.sum(score)
    
