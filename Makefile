CXX      ?= g++
CC       ?= gcc
CXXFLAGS  = -O2 -Wall -Wextra -std=c++17
CFLAGS    = -O2 -Wall
LDFLAGS   = -lrtlsdr -lpthread -lm

# macOS: Homebrew installs to /opt/homebrew (Apple Silicon) or /usr/local (Intel)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)
  CXXFLAGS += -I$(BREW_PREFIX)/include
  CFLAGS   += -I$(BREW_PREFIX)/include
  LDFLAGS  += -L$(BREW_PREFIX)/lib
endif

TARGET = lora_rx

SRCS_CXX = lora_rx.cpp
SRCS_C   = kiss_fft.c
OBJS     = $(SRCS_CXX:.cpp=.o) $(SRCS_C:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
