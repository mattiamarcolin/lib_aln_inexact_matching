TARGET_EXEC ?= example
CC = gcc
CFLAGS = -Wno-unused-result
BUILD_DIR ?= build
SRC_DIR ?= src

SRCS := $(shell find $(SRC_DIR) \( -name *.c \) )
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

OPTIM ?= -O3
LDXXFLAGS = -lpthread -lz 

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(OPTIM) -c $< -o $@

$(TARGET_EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OPTIM) $(OBJS) -o $@ $(LDXXFLAGS)

MKDIR_P ?= mkdir -p
