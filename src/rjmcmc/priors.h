// Copyright (C) 2012  Andrew H. Chan, Paul A. Jenkins, Yun S. Song
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// <http://www.gnu.org/licenses/>
//
// 

#ifndef LDHELMET_RJMCMC_PRIORS_H_
#define LDHELMET_RJMCMC_PRIORS_H_

#include <gsl/gsl_sf_gamma.h>

#include "rjmcmc/proposals.h"

double const PI = 3.14159265358979323846264338;

class AcceptanceRatio {
 public:
  virtual ~AcceptanceRatio() { }

  virtual double LogMHChange(double proposed_log_lk,
                             double cur_log_lk,
                             double new_rate,
                             double old_rate) const = 0;

  virtual double LogMHExtend(double proposed_log_lk,
                             double cur_log_lk) const = 0;

  virtual double LogMHSplit(double proposed_log_lk,
                            double cur_log_lk,
                            size_t num_snps,
                            size_t orig_num_change_points,
                            double new_left_rate,
                            double new_split_rate,
                            double old_rate,
                            Proposals const &proposals_) const = 0;

  virtual double LogMHMerge(double proposed_log_lk,
                            double cur_log_lk,
                            size_t num_snps,
                            size_t orig_num_change_points,
                            double left_rate,
                            double merge_rate,
                            double new_rate,
                            Proposals const &proposals_) const = 0;
};

class GammaPrior : public AcceptanceRatio {
 public:
  double block_penalty_;
  double alpha_, beta_;

  GammaPrior(double block_penalty, double alpha, double beta)
      : block_penalty_(block_penalty), alpha_(alpha), beta_(beta) { }

  double LogMHChange(double proposed_log_lk,
                     double cur_log_lk,
                     double new_rate,
                     double old_rate) const {
    double log_mh = proposed_log_lk - cur_log_lk
                  + alpha_ * (std::log(new_rate / old_rate))
                  - beta_ * (new_rate - old_rate);
    return log_mh;
  }

  double LogMHExtend(double proposed_log_lk, double cur_log_lk) const {
    double log_mh = proposed_log_lk - cur_log_lk;
    return log_mh;
  }

  double LogMHSplit(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double new_left_rate,
                    double new_split_rate,
                    double old_rate,
                    Proposals const &proposals_) const {
    double log_mh =
      proposed_log_lk - cur_log_lk
      - std::log(static_cast<double>((num_snps - 2) * (num_snps - 2))
                 / static_cast<double>((num_snps - orig_num_change_points)
                                       * (orig_num_change_points - 1)))
      - block_penalty_
      + alpha_ * std::log(beta_)
      - gsl_sf_lngamma(alpha_)
      + (alpha_ - 1.0)
          * std::log(new_left_rate * new_split_rate / old_rate)
      - beta_ * (new_left_rate + new_split_rate - old_rate)
      + std::log(proposals_.dist_.merge / proposals_.dist_.split)
      + 2.0 * std::log(new_left_rate + new_split_rate)
      - std::log(old_rate);
    return log_mh;
  }

  double LogMHMerge(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double left_rate,
                    double merge_rate,
                    double new_rate,
                    Proposals const &proposals_) const {
    double log_mh =
        proposed_log_lk - cur_log_lk
        - std::log(static_cast<double>((num_snps - orig_num_change_points)
                                       * (orig_num_change_points - 2))
                   / static_cast<double>((num_snps - 2) * (num_snps - 2)))
        + block_penalty_
        - alpha_ * std::log(beta_)
        + gsl_sf_lngamma(alpha_)
        - (alpha_ - 1.0) * std::log(left_rate * merge_rate / new_rate)
        + beta_ * (left_rate + merge_rate - new_rate)
        + std::log(proposals_.dist_.split / proposals_.dist_.merge)
        + std::log(new_rate)
        - 2.0 * std::log(left_rate + merge_rate);
    return log_mh;
  }
};

class LogNormalPrior : public AcceptanceRatio {
 public:
  double block_penalty_;
  double mu_, sigma_;

  LogNormalPrior(double block_penalty, double mu, double sigma)
      : block_penalty_(block_penalty), mu_(mu), sigma_(sigma) { }

  double LogMHChange(double proposed_log_lk,
                     double cur_log_lk,
                     double new_rate,
                     double old_rate) const {
    double log_mh =
        proposed_log_lk - cur_log_lk
        + ((std::log(old_rate) - mu_) * (std::log(old_rate) - mu_)
             - (std::log(new_rate) - mu_) * (std::log(new_rate) - mu_))
        / (2.0 * sigma_ * sigma_);
    return log_mh;
  }

  double LogMHExtend(double proposed_log_lk, double cur_log_lk) const {
    double log_mh = proposed_log_lk - cur_log_lk;
    return log_mh;
  }

  double LogMHSplit(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double new_left_rate,
                    double new_split_rate,
                    double old_rate,
                    Proposals const &proposals_) const {
    double log_mh =
      proposed_log_lk - cur_log_lk
      - std::log(static_cast<double>((num_snps - 2) * (num_snps - 2))
                   / static_cast<double>((num_snps - orig_num_change_points)
                                         * (orig_num_change_points - 1)))
      - block_penalty_
      + std::log(old_rate / (new_left_rate * new_split_rate))
      + std::log(std::pow(2.0 * PI * sigma_ * sigma_, 0.5))
      + (std::log(old_rate - mu_)*std::log(old_rate - mu_)
           - std::log(new_left_rate - mu_) * std::log(new_left_rate - mu_)
           - std::log(new_split_rate - mu_) * std::log(new_split_rate - mu_))
      / (2 * sigma_ * sigma_)
      + std::log(proposals_.dist_.merge / proposals_.dist_.split)
      + 2.0 * std::log(new_left_rate + new_split_rate)
      - std::log(old_rate);
    return log_mh;
  }

  double LogMHMerge(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double left_rate,
                    double merge_rate,
                    double new_rate,
                    Proposals const &proposals_) const {
    double log_mh =
      proposed_log_lk - cur_log_lk
      - std::log(static_cast<double>((num_snps - orig_num_change_points)
                                     * (orig_num_change_points - 2))
                 / static_cast<double>((num_snps - 2) * (num_snps - 2)))
      + block_penalty_
      - std::log(new_rate / (left_rate * merge_rate))
      - std::log(std::pow(2.0 * PI * sigma_ * sigma_, 0.5))
      - (std::log(new_rate - mu_) * std::log(new_rate - mu_)
           - std::log(left_rate - mu_) * std::log(left_rate - mu_)
           - std::log(merge_rate - mu_) * std::log(merge_rate - mu_))
      / (2 * sigma_ * sigma_)
      + std::log(proposals_.dist_.split / proposals_.dist_.merge)
      + std::log(new_rate)
      - 2.0 * std::log(left_rate + merge_rate);
    return log_mh;
  }
};

class ExponentialPrior : public AcceptanceRatio {
 public:
  double block_penalty_;
  double mean_;

  ExponentialPrior(double block_penalty, double mean)
      : block_penalty_(block_penalty), mean_(mean) { }

  double LogMHChange(double proposed_log_lk,
                     double cur_log_lk,
                     double new_rate,
                     double old_rate) const {
    double log_mh = proposed_log_lk - cur_log_lk
                  + std::log(new_rate / old_rate)
                  - (new_rate - old_rate) / mean_;
    return log_mh;
  }

  double LogMHExtend(double proposed_log_lk, double cur_log_lk) const {
    double log_mh = proposed_log_lk - cur_log_lk;
    return log_mh;
  }

    // Incorrect function from v1.7
  // double LogMHSplit(double proposed_log_lk,
  //                   double cur_log_lk,
  //                   size_t num_snps,
  //                   size_t orig_num_change_points,
  //                   double new_left_rate,
  //                   double new_split_rate,
  //                   double old_rate,
  //                   Proposals const &proposals_) const {
  //   double log_mh =
  //     proposed_log_lk - cur_log_lk
  //     - std::log(static_cast<double>((num_snps - 2) * (num_snps - 2))
  //                / static_cast<double>((num_snps - orig_num_change_points)
  //                                      * (orig_num_change_points - 1)))
  //     - block_penalty_
  //     - std::log(mean_)
  //     - (new_left_rate + new_split_rate - old_rate) / mean_
  //     + std::log(proposals_.dist_.merge / proposals_.dist_.split)
  //     + 2.0 * std::log(new_left_rate + new_split_rate)
  //     - std::log(old_rate);
  //   return log_mh;
  // }

  double LogMHSplit(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double new_left_rate,
                    double new_split_rate,
                    double old_rate,
                    Proposals const &proposals_) const {
    double log_mh =
      proposed_log_lk - cur_log_lk
      - std::log(static_cast<double>((orig_num_change_points + 1))
                 / static_cast<double>((num_snps - 2)))
      // - std::log(static_cast<double>((orig_num_change_points - 1) * (orig_num_change_points + 1))
      //            / static_cast<double>((num_snps - orig_num_change_points)
      //                                  * (num_snps - 2)))
      - block_penalty_
      - std::log(mean_)
      - (new_left_rate + new_split_rate - old_rate) / mean_
      + std::log(proposals_.dist_.merge / proposals_.dist_.split)
      + 2.0 * std::log(new_left_rate + new_split_rate)
      - std::log(old_rate);
    return log_mh;
  }
    
    // Incorrect function from v1.7
  // double LogMHMerge(double proposed_log_lk,
  //                   double cur_log_lk,
  //                   size_t num_snps,
  //                   size_t orig_num_change_points,
  //                   double left_rate,
  //                   double merge_rate,
  //                   double new_rate,
  //                   Proposals const &proposals_) const {
  //   double log_mh =
  //     proposed_log_lk - cur_log_lk
  //     - std::log(static_cast<double>((num_snps - orig_num_change_points)
  //                                    * (orig_num_change_points - 2))
  //                  / static_cast<double>((num_snps - 2) * (num_snps - 2)))
  //     + block_penalty_
  //     + std::log(mean_)
  //     + (left_rate + merge_rate - new_rate) / mean_
  //     + std::log(proposals_.dist_.split / proposals_.dist_.merge)
  //     + std::log(new_rate)
  //     - 2.0 * std::log(left_rate + merge_rate);
  //   return log_mh;
  // }

  double LogMHMerge(double proposed_log_lk,
                    double cur_log_lk,
                    size_t num_snps,
                    size_t orig_num_change_points,
                    double left_rate,
                    double merge_rate,
                    double new_rate,
                    Proposals const &proposals_) const {
    double log_mh =
      proposed_log_lk - cur_log_lk
      - std::log(static_cast<double>((num_snps - 2))
                   / static_cast<double>((orig_num_change_points)))
      // - std::log(static_cast<double>((num_snps - orig_num_change_points + 1)
      //                                * (num_snps - 2))
      //              / static_cast<double>((orig_num_change_points - 2) * (orig_num_change_points)))
      + block_penalty_
      + std::log(mean_)
      + (left_rate + merge_rate - new_rate) / mean_
      + std::log(proposals_.dist_.split / proposals_.dist_.merge)
      + std::log(new_rate)
      - 2.0 * std::log(left_rate + merge_rate);
    return log_mh;
  }
};

#endif
