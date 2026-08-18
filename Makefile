CXX      ?= g++
CC       ?= gcc
CXXFLAGS  = -O3 -fno-math-errno -Wall -Wextra -std=c++17
CFLAGS    = -O3 -fno-math-errno -Wall
LDFLAGS   = -lrtlsdr -lpthread -lm

# macOS: Homebrew installs to /opt/homebrew (Apple Silicon) or /usr/local (Intel)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)
  CXXFLAGS += -I$(BREW_PREFIX)/include
  CFLAGS   += -I$(BREW_PREFIX)/include
  LDFLAGS  += -L$(BREW_PREFIX)/lib
endif

# Optional HackRF input backend (-H): compiled in when libhackrf is present
HACKRF_LIBS := $(shell pkg-config --libs libhackrf 2>/dev/null)
ifneq ($(HACKRF_LIBS),)
  CXXFLAGS += -DHAVE_HACKRF $(shell pkg-config --cflags libhackrf 2>/dev/null)
  LDFLAGS  += $(HACKRF_LIBS)
endif

TARGET = lora_rx
TXGEN  = lora_tx_gen

SRCS_CXX = lora_rx.cpp
SRCS_C   = kiss_fft.c
OBJS     = $(SRCS_CXX:.cpp=.o) $(SRCS_C:.c=.o)

all: $(TARGET) $(TXGEN)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(TXGEN): lora_tx_gen.cpp lora_common.h
	$(CXX) $(CXXFLAGS) -o $@ lora_tx_gen.cpp -lm

lora_rx.o: lora_rx.cpp lora_common.h

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

check: $(TARGET) $(TXGEN)
	./run_tests.sh

clean:
	rm -f $(OBJS) $(TARGET) $(TXGEN)

.PHONY: all clean check
