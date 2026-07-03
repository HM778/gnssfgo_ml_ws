#ifndef LAMBDA_HPP
#define LAMBDA_HPP

#include <Eigen/Eigen>

/*
    * @brief LAMBDA integer least squares with decorrelation and search
    * @param n number of ambiguities
    * @param m number of candidates (expects m>=2 in current usage)
    * @param a float ambiguities (n x 1)
    * @param Q covariance (n x n), column-major
    * @param F integer candidates (n x m), column-major
    * @param S residuals (m x 1), S[0]=best, S[1]=second-best
    */
bool lambda(int n, int m, const double* a, const double* Q, double* F, double* S)
{
    if (n <= 0 || m <= 0 || !a || !Q || !F || !S)
    {
        return false;
    }

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> Qm(Q, n, n);
    Eigen::Map<const Eigen::VectorXd> a_vec(a, n);

    Eigen::LLT<Eigen::MatrixXd> llt(Qm);
    if (llt.info() != Eigen::Success)
    {
        return false;
    }

    Eigen::MatrixXd Lc = llt.matrixL();
    Eigen::VectorXd d = Lc.diagonal();
    if ((d.array() <= 0.0).any())
    {
        return false;
    }

    Eigen::VectorXd D = d.array().square();
    Eigen::MatrixXd L = Lc * d.cwiseInverse().asDiagonal();
    Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(n, n);

    int k = n - 2;
    while (k >= 0)
    {
        for (int j = k + 1; j < n; ++j)
        {
            double mu = std::round(L(j, k));
            if (mu != 0.0)
            {
                L.row(j).head(k + 1) -= mu * L.row(k).head(k + 1);
                Z.col(j) -= mu * Z.col(k);
            }
        }

        double delta = D(k) + L(k + 1, k) * L(k + 1, k) * D(k + 1);
        if (delta < D(k + 1))
        {
            double lam = D(k) / delta;
            double eta = D(k + 1) * L(k + 1, k) / delta;
            D(k) = lam * D(k + 1);
            D(k + 1) = delta;

            for (int i = 0; i < k; ++i)
            {
                std::swap(L(k, i), L(k + 1, i));
            }
            L(k + 1, k) = eta;
            for (int i = k + 2; i < n; ++i)
            {
                std::swap(L(i, k), L(i, k + 1));
            }
            Z.col(k).swap(Z.col(k + 1));
            k = std::min(k + 1, n - 2);
        }
        else
        {
            --k;
        }
    }

    Eigen::VectorXd z = Z.transpose() * a_vec;

    Eigen::MatrixXd Fz = Eigen::MatrixXd::Zero(n, m);
    Eigen::VectorXd Svec = Eigen::VectorXd::Constant(m, std::numeric_limits<double>::infinity());

    Eigen::VectorXd zcond = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd zhat = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd step = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dist = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dist_sum = Eigen::VectorXd::Zero(n);

    auto insert_candidate = [&](const Eigen::VectorXd& cand, double s) {
        if (s >= Svec(m - 1))
        {
            return;
        }
        int idx = m - 1;
        while (idx > 0 && s < Svec(idx - 1))
        {
            Svec(idx) = Svec(idx - 1);
            Fz.col(idx) = Fz.col(idx - 1);
            --idx;
        }
        Svec(idx) = s;
        Fz.col(idx) = cand;
    };

    int kk = n - 1;
    zcond(kk) = z(kk);
    zhat(kk) = std::round(zcond(kk));
    step(kk) = (zcond(kk) - zhat(kk) >= 0.0) ? 1.0 : -1.0;
    dist(kk) = (zcond(kk) - zhat(kk)) * (zcond(kk) - zhat(kk)) / D(kk);
    dist_sum(kk) = dist(kk);

    double maxdist = std::numeric_limits<double>::infinity();

    while (true)
    {
        double newdist = (kk == n - 1) ? dist(kk) : dist(kk) + dist_sum(kk + 1);
        if (newdist < maxdist)
        {
            if (kk == 0)
            {
                insert_candidate(zhat, newdist);
                maxdist = Svec(m - 1);

                zhat(0) += step(0);
                double y = zcond(0) - zhat(0);
                step(0) = (step(0) > 0.0) ? -step(0) - 1.0 : -step(0) + 1.0;
                dist(0) = y * y / D(0);
            }
            else
            {
                --kk;
                zcond(kk) = z(kk);
                for (int i = kk + 1; i < n; ++i)
                {
                    zcond(kk) += L(i, kk) * (zcond(i) - zhat(i));
                }
                zhat(kk) = std::round(zcond(kk));
                double y = zcond(kk) - zhat(kk);
                step(kk) = (y >= 0.0) ? 1.0 : -1.0;
                dist(kk) = y * y / D(kk);
            }
        }
        else
        {
            if (kk == n - 1)
            {
                break;
            }
            ++kk;
            zhat(kk) += step(kk);
            double y = zcond(kk) - zhat(kk);
            step(kk) = (step(kk) > 0.0) ? -step(kk) - 1.0 : -step(kk) + 1.0;
            dist(kk) = y * y / D(kk);
        }

        if (kk < n - 1)
        {
            dist_sum(kk) = dist(kk) + dist_sum(kk + 1);
        }
    }

    Eigen::MatrixXd ZinvT = Z.transpose().inverse();
    Eigen::MatrixXd F_a = ZinvT * Fz;

    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            F[i + n * j] = F_a(i, j);
        }
        S[j] = Svec(j);
    }

    return true;
}

#endif // LAMBDA_HPP