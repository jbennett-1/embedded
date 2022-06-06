CC := arm-none-eabi-gcc
LD := arm-none-eabi-ld
OBJCOPY := arm-none-eabi-objcopy
OBJDUMP := arm-none-eabi-objdump
SIZE := arm-none-eabi-size

CFLAGS := -O0 -fno-pie -ffreestanding -nostartfiles -fno-stack-protector -g3 -mthumb -Wall -mfloat-abi=hard -mfpu=fpv4-sp-d16 -march=armv7e-m
CFLAGS += -I./Include -I./Include/dsp -ffunction-sections -mno-unaligned-access
LDFLAGS := -Wl,--gc-sections

ODIR := obj
SDIR := src
SLIB := lib

OBJS = \
	startup_ARMCM4.o \
	system_ARMCM4.o \
	timeOptimization.o \
	eig_vec_decomp_micro.o \
	first.o

LIBS = \
	libCMSISDSPTransform.a \
	libCMSISDSPCommon.a \
	libCMSISDSPMatrix.a \
	libCMSISDSPComplexMath.a \
	libCMSISDSPFastMath.a \
	libCMSISDSPBasicMath.a \

LIB = $(patsubst %,$(SLIB)/%,$(LIBS))

OBJ = $(patsubst %,$(ODIR)/%,$(OBJS))

$(ODIR)/%.o: $(SLIB)/%.a
	$(LD) $(LDFLAGS) -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.s
	$(CC) $(CFLAGS) -c -g -o $@ $^ 

all: emb

emb: $(OBJ)
	$(CC) $(CFLAGS) obj/* $(LIB) -Tgcc_arm.ld $(LDFLAGS) -o embedded.img -lm
	cp embedded.img embedded.elf
	$(OBJCOPY) -O binary embedded.img
	$(SIZE) embedded.elf

clean:
	rm -f obj/*
	rm -f embedded.elf

debug:
	screen -S openocd -d -m openocd -f /home/juliabennett/Desktop/openocd/openocd-code/tcl/board/microchip_same54_xplained_pro.cfg
	TERM=xterm gdb-multiarch -x gdb_init_prot_mode.txt embedded.elf

run:
	screen -S openocd -d -m openocd -f /home/juliabennett/Desktop/openocd/openocd-code/tcl/board/microchip_same54_xplained_pro.cfg -kernel embedded.elf -serial null -monitor stdio

disassemble:
	$(OBJDUMP) -d embedded.elf


