SHELL = /bin/sh

ifeq (${OS},Windows_NT)
  PLATFORM = win32
else
  PLATFORM = $(shell uname|tr '[:upper:]' '[:lower:]')
endif
include ../Makefile.$(PLATFORM)
