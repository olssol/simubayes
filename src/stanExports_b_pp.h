// Generated by rstantools.  Do not edit by hand.

/*
    test is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    test is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with test.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_b_pp_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 14> locations_array__ =
  {" (found before start of program)",
  " (in 'b_pp', line 10, column 6 to column 35)",
  " (in 'b_pp', line 13, column 6 to column 64)",
  " (in 'b_pp', line 17, column 8 to column 42)",
  " (in 'b_pp', line 18, column 8 to column 29)",
  " (in 'b_pp', line 2, column 6 to column 12)",
  " (in 'b_pp', line 3, column 13 to column 14)",
  " (in 'b_pp', line 3, column 6 to column 19)",
  " (in 'b_pp', line 4, column 13 to column 14)",
  " (in 'b_pp', line 4, column 6 to column 19)",
  " (in 'b_pp', line 5, column 6 to column 23)",
  " (in 'b_pp', line 6, column 6 to column 23)",
  " (in 'b_pp', line 7, column 6 to column 25)",
  " (in 'b_pp', line 13, column 13 to column 14)"};
#include <stan_meta_header.hpp>
class model_b_pp final : public model_base_crtp<model_b_pp> {
private:
  int k;
  Eigen::Matrix<double,-1,1> N0_data__;
  Eigen::Matrix<double,-1,1> y0_data__;
  double c;
  double d;
  double a_0;
  Eigen::Map<Eigen::Matrix<double,-1,1>> N0{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,1>> y0{nullptr, 0};
public:
  ~model_b_pp() {}
  model_b_pp(stan::io::var_context& context__, unsigned int
             random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_b_pp_namespace::model_b_pp";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 5;
      context__.validate_dims("data initialization", "k", "int",
        std::vector<size_t>{});
      k = std::numeric_limits<int>::min();
      current_statement__ = 5;
      k = context__.vals_i("k")[(1 - 1)];
      current_statement__ = 6;
      stan::math::validate_non_negative_index("N0", "k", k);
      current_statement__ = 7;
      context__.validate_dims("data initialization", "N0", "double",
        std::vector<size_t>{static_cast<size_t>(k)});
      N0_data__ = Eigen::Matrix<double,-1,1>::Constant(k,
                    std::numeric_limits<double>::quiet_NaN());
      new (&N0) Eigen::Map<Eigen::Matrix<double,-1,1>>(N0_data__.data(), k);
      {
        std::vector<local_scalar_t__> N0_flat__;
        current_statement__ = 7;
        N0_flat__ = context__.vals_r("N0");
        current_statement__ = 7;
        pos__ = 1;
        current_statement__ = 7;
        for (int sym1__ = 1; sym1__ <= k; ++sym1__) {
          current_statement__ = 7;
          stan::model::assign(N0, N0_flat__[(pos__ - 1)],
            "assigning variable N0", stan::model::index_uni(sym1__));
          current_statement__ = 7;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 8;
      stan::math::validate_non_negative_index("y0", "k", k);
      current_statement__ = 9;
      context__.validate_dims("data initialization", "y0", "double",
        std::vector<size_t>{static_cast<size_t>(k)});
      y0_data__ = Eigen::Matrix<double,-1,1>::Constant(k,
                    std::numeric_limits<double>::quiet_NaN());
      new (&y0) Eigen::Map<Eigen::Matrix<double,-1,1>>(y0_data__.data(), k);
      {
        std::vector<local_scalar_t__> y0_flat__;
        current_statement__ = 9;
        y0_flat__ = context__.vals_r("y0");
        current_statement__ = 9;
        pos__ = 1;
        current_statement__ = 9;
        for (int sym1__ = 1; sym1__ <= k; ++sym1__) {
          current_statement__ = 9;
          stan::model::assign(y0, y0_flat__[(pos__ - 1)],
            "assigning variable y0", stan::model::index_uni(sym1__));
          current_statement__ = 9;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 10;
      context__.validate_dims("data initialization", "c", "double",
        std::vector<size_t>{});
      c = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 10;
      c = context__.vals_r("c")[(1 - 1)];
      current_statement__ = 10;
      stan::math::check_greater_or_equal(function__, "c", c, 0);
      current_statement__ = 11;
      context__.validate_dims("data initialization", "d", "double",
        std::vector<size_t>{});
      d = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 11;
      d = context__.vals_r("d")[(1 - 1)];
      current_statement__ = 11;
      stan::math::check_greater_or_equal(function__, "d", d, 0);
      current_statement__ = 12;
      context__.validate_dims("data initialization", "a_0", "double",
        std::vector<size_t>{});
      a_0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 12;
      a_0 = context__.vals_r("a_0")[(1 - 1)];
      current_statement__ = 12;
      stan::math::check_greater_or_equal(function__, "a_0", a_0, 0);
      current_statement__ = 13;
      stan::math::validate_non_negative_index("logL", "k", k);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1;
  }
  inline std::string model_name() const final {
    return "model_b_pp";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_b_pp_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      local_scalar_t__ theta = DUMMY_VAR__;
      current_statement__ = 1;
      theta = in__.template read_constrain_lub<local_scalar_t__,
                jacobian__>(0, 1, lp__);
      Eigen::Matrix<local_scalar_t__,-1,1> logL =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(k, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(logL,
        stan::math::add(stan::math::multiply(y0, stan::math::log(theta)),
          stan::math::multiply(stan::math::subtract(N0, y0),
            stan::math::log1m(theta))), "assigning variable logL");
      {
        current_statement__ = 3;
        lp_accum__.add(stan::math::beta_lpdf<false>(theta, c, d));
        current_statement__ = 4;
        lp_accum__.add(stan::math::multiply(a_0, logL));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_b_pp_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      double theta = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      theta = in__.template read_constrain_lub<local_scalar_t__,
                jacobian__>(0, 1, lp__);
      Eigen::Matrix<double,-1,1> logL =
        Eigen::Matrix<double,-1,1>::Constant(k,
          std::numeric_limits<double>::quiet_NaN());
      out__.write(theta);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 2;
      stan::model::assign(logL,
        stan::math::add(stan::math::multiply(y0, stan::math::log(theta)),
          stan::math::multiply(stan::math::subtract(N0, y0),
            stan::math::log1m(theta))), "assigning variable logL");
      if (emit_transformed_parameters__) {
        out__.write(logL);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ theta = DUMMY_VAR__;
      current_statement__ = 1;
      theta = in__.read<local_scalar_t__>();
      out__.write_free_lub(0, 1, theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "theta", "double",
        std::vector<size_t>{});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ theta = DUMMY_VAR__;
      current_statement__ = 1;
      theta = context__.vals_r("theta")[(1 - 1)];
      out__.write_free_lub(0, 1, theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"theta"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"logL"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(k)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "theta");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= k; ++sym1__) {
        param_names__.emplace_back(std::string() + "logL" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "theta");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= k; ++sym1__) {
        param_names__.emplace_back(std::string() + "logL" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"logL\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(k) + "},\"block\":\"transformed_parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"logL\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(k) + "},\"block\":\"transformed_parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = 1;
    const size_t num_transformed = emit_transformed_parameters * (k);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = 1;
    const size_t num_transformed = emit_transformed_parameters * (k);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_b_pp_namespace::model_b_pp;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_b_pp_namespace::profiles__;
}
#endif
#endif
