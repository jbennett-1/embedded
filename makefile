CC := arm-none-eabi-gcc
LD := arm-none-eabi-ld
OBJCOPY := arm-none-eabi-objcopy
OBJDUMP := arm-none-eabi-objdump
SIZE := arm-none-eabi-size

CFLAGS := -O0 -ffreestanding -fno-pie -fno-stack-protector -g3 -march=armv7e-m -mthumb -Wall -mfloat-abi=hard -mfpu=fpv4-sp-d16
CFLAGS += -I/home/juliabennett/Desktop/embedded/Include
LDFLAGS := -L/home/juliabennett/Desktop/embedded/lib 

ODIR := obj
SDIR := src
SLIB := lib

OBJS = \
	startup_ARMCM7.o \
	system_ARMCM7.o \
	mainArm.o	\
	timeOptimization.o \
	eig_vec_decomp_micro.o

LIBS = \
	libCMSISDSPTransform.a \
	libCMSISDSPBasicMath.a \
	libCMSISDSPMatrix.a \
	libCMSISDSPCommon.a \

LIB = $(patsubst %,$(SLIB)/%,$(LIBS))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJS))

$(ODIR)/%.o: $(SLIB)/%.a
	$(LD) $(LDFLAGS) $(LDLIBS) -o $@ $^ 

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.s
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.S
	$(CC) $(CFLAGS) -c -g -o $@ $^

all: emb

emb: $(OBJ) 
	$(LD) obj/* -Tgcc_arm.ld $(LIB) -o embedded.img 
	cp embedded.img embedded.elf
	$(OBJCOPY) -O binary embedded.img
	$(SIZE) embedded.elf

clean:
	rm -f obj/*
	rm -f embedded.elf

debug:  
	screen -S openocd -d -m openocd -f /home/juliabennett/Desktop/openocd/tcl/target/atsame5x.cfg
	gdb-multiarch embedded.elf 
	TERM=xterm gdb-multiarch -x gdb_init_prot_mode.txt

run:
	screen -S openocd -d -m openocd -f /home/juliabennett/Desktop/openocd/tcl/target/atsame5x.cfg -kernel embedded.elf -serial null -monitor stdio

disassemble:
	$(OBJDUMP) -d embedded.elf


