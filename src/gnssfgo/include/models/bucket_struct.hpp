#ifndef BUCKET_STRUCT_HPP
#define BUCKET_STRUCT_HPP
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <deque>
#include "../datatype.h"


// 每个桶统计TR因子的残差分布
struct BucketStat
{
    int count = 0;
    int success_count = 0;

    // 更新一个新残差
    void update(double residual)
    {
        count++;
        // printf("new update sample: %f\n",residual);
        if (residual < 0.005)
        {
            success_count++;
        }
    }

    void remove(double residual)
    {
        if (count > 0)
        {
            count--;
        }
        if (residual < 0.05 && success_count > 0)
        {
            success_count--;
        }
    }

    bool empty() const
    {
        return count <= 0;
    }

};

struct TimeDiscountSample
{
    int epoch_index = -1;
    double gap_epoch = 0.0;
    double residual = 0.0;
    double delta_pos = 0.0;
};

// 独立时间间隔处理结构：不做时间分桶，只聚合所有时间样本学习折扣函数
struct TimeDiscountStat
{
    std::deque<TimeDiscountSample> samples;
    double a,b;
    double mean_delta_pos=0.0;
    double discount = 1.0;

    Eigen::Matrix<double, Eigen::Dynamic, 2> toObsMatrix() const{
        int rows = samples.size();
        Eigen::Matrix<double, Eigen::Dynamic, 2> mat(rows, 2);
    
        int i = 0;
        for (const auto& item : samples) {
            mat(i, 0) = item.gap_epoch;
            mat(i, 1) = std::pow(item.delta_pos, 2);
            ++i;
        }
        return mat;
    }
    Eigen::Matrix<double, Eigen::Dynamic, 1> toResMatrix() const{
        int rows = samples.size();
        Eigen::Matrix<double, Eigen::Dynamic, 1> mat(rows, 1);
    
        int i = 0;
        for (const auto& sample : samples) {
            mat(i, 0) = sample.residual;
            i++;
        }
        return mat;
    }

    void rebuildModel()
    {
        if (samples.size() < 2)
        {
            a = 0.0;
            b = 1.0;
            mean_delta_pos = 0.0;
            discount = 1.0;
            return;
        }

        auto obs_mat = toObsMatrix();
        auto r_vec = toResMatrix(); 
        Eigen::Matrix<double, 2, 1> coeffs = (obs_mat.transpose() * obs_mat).ldlt().solve(obs_mat.transpose() * r_vec);
        a = std::max(0.0, coeffs(0));  // 斜率
        b = coeffs(1);  // 截距
    }

    void update(int epoch_index, double gap_epoch, double residual, double delta_pos)
    {
        samples.push_back({epoch_index, gap_epoch, residual, delta_pos});
        rebuildModel();
    }

    void update(double gap_epoch,double residual, double delta_pos)
    {
        update(-1, gap_epoch, residual, delta_pos);
    }

    void removeEpoch(int epoch_index)
    {
        if (epoch_index < 0)
        {
            return;
        }

        const auto old_size = samples.size();
        samples.erase(std::remove_if(samples.begin(), samples.end(),
                                     [epoch_index](const TimeDiscountSample& sample)
                                     {
                                         return sample.epoch_index == epoch_index;
                                     }),
                      samples.end());

        if (samples.size() != old_size)
        {
            rebuildModel();
        }
    }

    void clear()
    {
        samples.clear();
        a = 0.0;
        b = 1.0;
        mean_delta_pos = 0.0;
        discount = 1.0;
    }

    double getDiscount(double gap_epoch,double delta_pos) const
    {
        if (!std::isfinite(gap_epoch) || gap_epoch < 0.0)
        {
            return 1.0; // 对于无效的时间间隔，返回默认折扣
        }
        double discount = a * gap_epoch + b * pow(delta_pos, 2); // 线性模型计算折扣
        // printf("recent learning | a: %.3f, b: %.3f, gap_epoch: %.3f, delta_pos: %.3f, raw_discount: %.3f\n", a, b, gap_epoch, delta_pos, discount);
        if (discount < 0.0)
        {
            discount = 0.001; // 设置折扣下限，避免过度惩罚
        }

        discount = discount * 100.0; // [0,1]

        if (discount >= 1.0)
        {
            discount = 0.99;
        }
        // 返回真正的折扣因子，范围限定在(0, 1]，避免把长基线/大位移样本反向放大。
        return (1.0-discount);
    }
    
};
// 桶键，根据TR因子的特征计算得到
struct BucketKey
{
    int snr_min_bin;   
    int snr_max_bin;   
    int elev_min_bin;    
    int elev_max_bin;    
    int delta_elev_bin;
    int delta_azm_bin;      // 时间间隔单独处理用于计算权重

    bool operator==(const BucketKey& other) const
    {
        return snr_min_bin == other.snr_min_bin &&
               snr_max_bin == other.snr_max_bin &&
               elev_min_bin == other.elev_min_bin &&
               elev_max_bin == other.elev_max_bin &&
               delta_elev_bin == other.delta_elev_bin &&
               delta_azm_bin == other.delta_azm_bin;
    }
};


struct BucketKeyHash  
{
    std::size_t operator()(const BucketKey& k) const
    {
        // include all dimensions (including geo_elev_bin) to reduce hash collisions
        std::size_t h = static_cast<std::size_t>(k.snr_min_bin);
        h = h * 131 + static_cast<std::size_t>(k.snr_max_bin);
        h = h * 131 + static_cast<std::size_t>(k.elev_min_bin);
        h = h * 131 + static_cast<std::size_t>(k.elev_max_bin);
        h = h * 131 + static_cast<std::size_t>(k.delta_elev_bin);
        h = h * 131 + static_cast<std::size_t>(k.delta_azm_bin);
        return h;
    }
};

class BucketManager
{
public:
    std::unordered_map<int, BucketStat> snr_min_buckets;
    std::unordered_map<int, BucketStat> snr_max_buckets;
    std::unordered_map<int, BucketStat> elev_min_buckets;
    std::unordered_map<int, BucketStat> elev_max_buckets;
    std::unordered_map<int, BucketStat> delta_elev_buckets;
    std::unordered_map<int, BucketStat> delta_azm_buckets;
    TimeDiscountStat time_discount_stat;

    // 更新一个TR因子的残差
    void update(const BucketKey& key, double residual)
    {
        update(key, residual, residual < 0.05);
    }

    void update(const BucketKey& key, double residual, bool success)
    {
        snr_min_buckets[key.snr_min_bin].update(residual);
        snr_max_buckets[key.snr_max_bin].update(residual);
        elev_min_buckets[key.elev_min_bin].update(residual);
        elev_max_buckets[key.elev_max_bin].update(residual);
        delta_elev_buckets[key.delta_elev_bin].update(residual);
        delta_azm_buckets[key.delta_azm_bin].update(residual);
        (void)success;
    }

    void remove(const BucketKey& key, double residual)
    {
        remove(key, residual, residual < 0.05);
    }

    void remove(const BucketKey& key, double residual, bool success)
    {
        removeFromBucket(snr_min_buckets, key.snr_min_bin, residual, success);
        removeFromBucket(snr_max_buckets, key.snr_max_bin, residual, success);
        removeFromBucket(elev_min_buckets, key.elev_min_bin, residual, success);
        removeFromBucket(elev_max_buckets, key.elev_max_bin, residual, success);
        removeFromBucket(delta_elev_buckets, key.delta_elev_bin, residual, success);
        removeFromBucket(delta_azm_buckets, key.delta_azm_bin, residual, success);
    }

    void updateTimeDiscount(int epoch_index, double gap_epoch, double residual,double delta_pos)
    {
        time_discount_stat.update(epoch_index, gap_epoch, residual, delta_pos);
    }

    void updateTimeDiscount(double gap_epoch, double residual,double delta_pos)
    {
        time_discount_stat.update(gap_epoch, residual, delta_pos);
    }

    void removeTimeDiscountEpoch(int epoch_index)
    {
        time_discount_stat.removeEpoch(epoch_index);
    }

    void clear()
    {
        snr_min_buckets.clear();
        snr_max_buckets.clear();
        elev_min_buckets.clear();
        elev_max_buckets.clear();
        delta_elev_buckets.clear();
        delta_azm_buckets.clear();
        time_discount_stat.clear();
    }

private:
    static void removeFromBucket(std::unordered_map<int, BucketStat>& bucket_map, int bin, double residual, bool success)
    {
        auto it = bucket_map.find(bin);
        if (it == bucket_map.end())
        {
            return;
        }

        it->second.remove(residual);
        if (it->second.empty())
        {
            bucket_map.erase(it);
        }
        (void)success;
    }

    static double getBinReliability(const std::unordered_map<int, BucketStat>& bucket_map, int bin, double default_prob = 0.5)
    {
        const auto it = bucket_map.find(bin);
        if (it == bucket_map.end())
        {
            return default_prob;
        }

        const BucketStat& stat = it->second;
        if (stat.count <= 0)
        {
            return default_prob;
        }

        const double probability = static_cast<double>(stat.success_count) / static_cast<double>(stat.count);
        return std::max(0.0, std::min(1.0, probability));
    }

    double getCompositeReliability(const BucketKey& key) const
    {
        const double p_snr_min = getBinReliability(snr_min_buckets, key.snr_min_bin);
        const double p_snr_max = getBinReliability(snr_max_buckets, key.snr_max_bin);
        const double p_elev_min = getBinReliability(elev_min_buckets, key.elev_min_bin);
        const double p_elev_max = getBinReliability(elev_max_buckets, key.elev_max_bin);
        const double p_delta_elev = getBinReliability(delta_elev_buckets, key.delta_elev_bin);
        const double p_delta_azm = getBinReliability(delta_azm_buckets, key.delta_azm_bin);

        // printf("Bucket reliabilities | SNR_min_bin %d: %.3f, SNR_max_bin %d: %.3f, Elev_min_bin %d: %.3f, Elev_max_bin %d: %.3f, Delta_Elev_bin %d: %.3f, Delta_Azm_bin %d: %.3f\n",
        //        key.snr_min_bin, p_snr_min,
        //        key.snr_max_bin, p_snr_max,
        //        key.elev_min_bin, p_elev_min,
        //        key.elev_max_bin, p_elev_max,
        //        key.delta_elev_bin, p_delta_elev,
        //        key.delta_azm_bin, p_delta_azm);
        const double probability = p_snr_min * p_snr_max * p_elev_min * p_elev_max * p_delta_elev * p_delta_azm;
        return std::max(0.0, std::min(1.0, probability));
    }

public:
    // 获取“整个键”的成功概率：各个分桶成功概率相乘
    double getReliability(const BucketKey& key) const
    {
        return getCompositeReliability(key);
    }

    // 为兼容旧调用，保留这个接口；这里也返回同一个成功概率。
    double getVariance(const BucketKey& key) const
    {
        return getCompositeReliability(key);
    }

    double getTimeDiscount(double gap_epoch,double delta_pos) const
    {
        return time_discount_stat.getDiscount(gap_epoch,delta_pos);
    }
};

// 桶索引器，根据TR因子的特征计算桶键
class BucketIndexer
{
public:
    int SnrToBin(double snr)
    {
        if (snr < 35) return 0;
        else if (snr < 40) return 1;
        else if (snr < 45) return 2;
        else return 3;
    }

    int EleToBin(double ele)
    {
        if (ele <= 10) return 0;
        else if (ele <= 30) return 1;
        else if (ele <= 60) return 2;
        else if (ele <= 90) return 3;
    }

    int DeltaEleToBin(double delta_ele)
    {
        if (delta_ele <= 10) return 0;
        else if (delta_ele <= 45) return 1;
        else if (delta_ele <= 80) return 2;
        else return 3;
    }
    
    int DeltaAzmToBin(double azm)
    {
        if (azm < 0.5) return 0;
        else if (azm < 1.0) return 1;
        else if (azm < 1.5) return 2;
        else return 3;
    }

    BucketKey makeKey(double min_snr, double max_snr, double min_ele, double max_ele, double azm)
    {
        BucketKey key;
        key.snr_min_bin = SnrToBin(min_snr);
        key.snr_max_bin = SnrToBin(max_snr);
        key.elev_min_bin = EleToBin(min_ele);
        key.elev_max_bin = EleToBin(max_ele);
        key.delta_elev_bin = DeltaEleToBin(max_ele - min_ele);
        key.delta_azm_bin = DeltaAzmToBin(azm);
        return key;
    }

};

#endif // BUCKET_STRUCT_HPP