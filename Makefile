# -------------- #
# -- Makefile -- #
# -------------- #

# Copyright (c) 2000-2013 Lionel Lacassagne
# with a little help of Stephane Piskorski and Joel Falcou


UNAME_S := $(shell uname -s)



	
# -- Lile list ----------
FILE =  img.c mouvement.c nrutil.c vnrutil.c mutil.c morpho.c\
		img_SIMD.c mouvement_SIMD.c morpho_SIMD.c test_morpho.c test_mouvement.c\
		util.c mynrutil.c morpho_optim.c morpho_pack_optim.c benchmark.c

ifneq ($(UNAME_S),Darwin)
	FILE += morpho_optim_omp.c morpho_main.c
else
	FILE += main.c
endif
 
# -- Paths ----------
SRC_PATH = src
OBJ_PATH = obj
EXE_PATH = exe
INC_PATH = include

# -- OS ----------
#OS = MACH_OSX

# -- Config ----------
# if CONFIG = CLI  (Command Line Interface, no Apple Framework)
CONFIG = CLI

# -- Macros ----------
CC = gcc
AR = ar -rc

# -- Flags ----------
C_DEBUG_FLAGS = -O3 -g
C_CC_FLAGS = -std=c99
C_OPTIMISATION_FLAGS = -O3 -fstrict-aliasing
C_ARCH_FLAGS = -mssse3 
ifneq ($(UNAME_S),Darwin)
	C_ARCH_FLAGS += -fopenmp
endif
#C_OS_FLAGS = -D$(OS)
C_CONFIG_FLAGS = -D$(CONFIG)
C_INC_FLAGS = -I$(INC_PATH)

CFLAGS =  $(C_CC_FLAGS) $(C_DEBUG_FLAGS)        $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS)
#CFLAGS = $(C_CC_FLAGS) $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS)
LDFLAGS = $(C_CC_FLAGS) $(C_OPTIMISATION_FLAGS) $(C_ARCH_FLAGS) $(C_OS_FLAGS) $(C_CONFIG_FLAGS) $(C_INC_FLAGS)

# -- Final product ----------
PRODUCT   = project

# -- src and obj List ----------
SRC = $(addprefix ${SRC_PATH}/, $(FILE))
OBJ = $(addprefix ${OBJ_PATH}/, $(addsuffix .o, $(basename $(FILE))))

# -- Base rules ----------
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

#-----Main rule ----------
$(EXE_PATH)/$(PRODUCT): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(OPTFLAGS) $(CFG) $(INC) $(LIB) -lm

# -- Other stuff ----------
depend:
	makedepend $(CFLAGS) -Y $(SRC)

clean:
	rm -f $(OBJ)
	rm -f ${EXE_PATH}/${PRODUCT}