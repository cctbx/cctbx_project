#include "iotbx/cppmtz.h"

namespace bpmtz = iotbx::mtz;

af::shared<std::string> bpmtz::Foo::value(){
   af::shared<std::string> answer;
   answer.push_back("TEST");
   return answer;
}

