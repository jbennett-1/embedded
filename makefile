CC := arm-none-eabi-gcc
LD := arm-none-eabi-ld
OBJCOPY := arm-none-eabi-objcopy
OBJDUMP := arm-none-eabi-objdump
SIZE := arm-none-eabi-size

CFLAGS := -O0 -ffreestanding -fno-pie -fno-stack-protector -g3 -march=armv7e-m -mthumb -Wall -mfloat-abi=hard -mfpu=fpv4-sp-d16
CFLAGS += -I/home/juliabennett/Desktop/embedded/Include -nostartfiles -ffunction-sections -mno-unaligned-access
LDFLAGS := -L/home/juliabennett/Desktop/embedded/lib -gc-sections

ODIR := obj
SDIR := src
SLIB := lib

OBJS = \
	startup_ARMCM4.o \
	system_ARMCM4.o \
	timeOptimization.o \
	eig_vec_decomp_micro.o \
	mainArm.o

LIBS = \
	libCMSISDSPTransform.a \
	libCMSISDSPCommon.a \
	libCMSISDSPMatrix.a \
	libCMSISDSPBasicMath.a \
	libCMSISDSPComplexMath.a

LIB = $(patsubst %,$(SLIB)/%,$(LIBS))

OBJ = $(patsubst %,$(ODIR)/%,$(OBJS))

$(ODIR)/%.o: $(SLIB)/%.a
	$(LD) $(LDFLAGS) $(LDLIBS) -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.s
	$(CC) $(CFLAGS) -c -g -o $@ $^

#$(ODIR)/%.o: $(SDIR)/%.S
#	$(CC) $(CFLAGS) -c -g -o $@ $^

all: emb

emb: $(OBJ)
	$(LD) obj/* -Tgcc_arm.ld $(LIB) -gc-sections -o embedded.img
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


