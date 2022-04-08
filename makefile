CC := arm-none-eabi-gcc
LD := arm-none-eabi-ld
OBJCOPY := arm-none-eabi-objcopy
OBJDUMP := arm-none-eabi-objdump
SIZE := arm-none-eabi-size

CFLAGS := -O0 -ffreestanding -fno-pie -fno-stack-protector -g3 -march=armv7-m -mthumb -Wall -mgeneral-regs-only 
CFLAGS += -I/home/juliabennett/Desktop/embeddedPCA/Include

ODIR := obj
SDIR := src
Lib := lib

Lib = \
      libCMSISDSPBasicMath.a \
      libCMSISDSPMatrix.a \
      libCMSISDSPTransform.a

OBJS = \
	main.o \
	startup_ARMCM7.o \
	system_ARMCM7.o \
	mainArm.o \
	timeOptimization.o

OBJ = $(patsubst %,$(ODIR)/%,$(OBJS))


$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.s
	$(CC) $(CFLAGS) -c -g -o $@ $^

$(ODIR)/%.o: $(SDIR)/%.S
	$(CC) $(CFLAGS) -c -g -o $@ $^

all: emb

emb: $(OBJ)
	$(LD) obj/* -Tgcc_arm.ld -o embedded.img
	cp embedded.img embedded.elf
	$(OBJCOPY) -O binary embedded.img
	$(SIZE) embedded.elf

clean:
	rm -f obj/*
	rm -f embedded.elf

debug:
	screen -S openocd -d -m openocd -f /usr/share/openocd/scripts/board/atmel_same54_explained.cfg 
	gdb-multiarch embedded.elf 
	TERM=xterm gdb-multiarch -x gdb_init_prot_mode.txt

run:
	screen -S openocd -d -m openocd -f /usr/share/openocd/scripts/board/atmel_same54_explained.cfg -kernel embedded.elf -serial null -monitor stdio

disassemble:
	$(OBJDUMP) -d embedded.elf


