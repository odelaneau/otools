#ifndef LIFTOVER_CHAINFILE_H
#define LIFTOVER_CHAINFILE_H

#include <utils/otools.h>

#include <containers/chain.h>
#include <containers/target.h>

namespace liftover {

std::map<std::string, Target> open_chainfile(std::string path);

}

#endif
