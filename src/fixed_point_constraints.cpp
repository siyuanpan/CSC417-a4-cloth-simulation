#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    P.resize(q_size - indices.size() * 3, q_size);
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(q_size - indices.size() * 3);

    unsigned int cnt = 0;
    for (unsigned int i = 0; i < q_size; i += 3)
    {
        if (std::find(indices.begin(), indices.end(), (unsigned int)(i / 3)) != indices.end())
            continue;
        triplets.push_back(T(cnt, i, 1.));
        triplets.push_back(T(cnt + 1, i + 1, 1.));
        triplets.push_back(T(cnt + 2, i + 2, 1.));
        cnt += 3;
    }
    P.setFromTriplets(triplets.begin(), triplets.end());
    

}
