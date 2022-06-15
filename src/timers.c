#define ARMCM4_FP
#include "ARMCM4_FP.h"
#include "core_cm4.h"

uint32_t ms_ticks = 0;

void SysTick_Handler(void)
{
    ms_ticks++;
}

int main(void)
{
    uint32_t returnCode;
    returnCode = SysTick_Config(SystemCoreClock / 100); //every 10 seconds
    
    if(returnCode != 0){
	SysTick_Handler();
    }
    while(1);
}
