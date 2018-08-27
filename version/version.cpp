#include <version/git-revision-internal.h>
#include <version/version.h>

#define STR_EXPAND(x) #x
#define STR(x) STR_EXPAND(x)
#define STRING_REVISION STR(REVISION)

namespace version {
const std::string revision = STRING_REVISION;
}
