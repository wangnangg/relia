#include <gtest/gtest.h>
#include "easylogging++.h"

int main(int ac, char* av[])
{
	el::Configurations conf("./test_logger.conf");
	el::Loggers::reconfigureAllLoggers(conf);
	testing::InitGoogleTest(&ac, av);
	return RUN_ALL_TESTS();
}