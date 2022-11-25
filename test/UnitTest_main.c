#include "UnitTest_LinearC.h"

int main()
{
    UnitTest_LinearC* unit_test;
    
    unit_test = UnitTest_LinearC_New();
    UnitTest_LinearC_TestAll(unit_test);
    UnitTest_LinearC_Del(unit_test);
    
    return 0;
}