// The BOUT++ half of the in-process tests.
// Compiles the unmodified elm-pb physics model.
#include <bout/bout.hxx>
#include <bout/monitor.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/solver.hxx>
#undef BOUTMAIN
#define BOUTMAIN(x) /* the example's main() is not wanted here */
#include ELMPB_SOURCE

#include "bout_connect/bout_runner.hpp"
#include "bout_connect/bout_interface.hpp"

namespace bout_connect
{
namespace
{
class EveryCall : public Monitor
{
public:
  explicit EveryCall(const std::function<void(double, int, int)>& f) : f_(f) {}
  int call(Solver*, BoutReal t, int iter, int nout) override
  {
    f_(t, iter, nout);
    return 0;
  }

private:
  const std::function<void(double, int, int)>& f_;
};
} // namespace

int run_elmpb(int argc, char** argv, const std::function<void()>& after_init,
              const SnapshotWork& work, bool& ran)
{
  const int rc = BoutInitialise(argc, argv);
  if (rc != 0)
    return rc;
  after_init();
  const int status = run_paused<ELMpb>(
    [&](Solver* solver, BoutReal t, int iter, int nout) {
      Captured c = capture(solver);
      c.data.xz_bilinear =
        [field = c.field](const std::vector<std::array<double, 3>>& q) {
          return xz_bilinear(field, q);
        };
      work(c.data, t, iter, nout);
    },
    ran);
  BoutFinalise();
  return status;
}

int run_elmpb_monitor(int argc, char** argv,
                      const std::function<void(double, int, int)>& on_call)
{
  const int rc = BoutInitialise(argc, argv);
  if (rc != 0)
    return rc;
  EveryCall mon(on_call);
  const int status = run_model<ELMpb>(mon);
  BoutFinalise();
  return status;
}
int option_int(const char* section, const char* key, const char* doc, int dflt)
{
  return Options::root()[section][key].doc(doc).withDefault(dflt);
}
double option_double(const char* section, const char* key, const char* doc,
                     double dflt)
{
  return Options::root()[section][key].doc(doc).withDefault(dflt);
}
bool option_bool(const char* section, const char* key, const char* doc,
                 bool dflt)
{
  return Options::root()[section][key].doc(doc).withDefault(dflt);
}

} // namespace bout_connect
