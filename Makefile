CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
OPENMP ?= -fopenmp
PREFIX ?= /usr/local
PYTHON ?= python3
SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp,build/%.o,$(SOURCES))

.PHONY: all test install clean
all: build/altann-core

build/%.o: src/%.cpp src/altann.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -MMD -MP -c $< -o $@

build/altann-core: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OPENMP) $(OBJECTS) $(LDFLAGS) -o $@

test: all
	$(PYTHON) -m unittest discover -s tests -v

# Keep the Python frontend and compatibility preprocessor together. The
# launcher resolves this directory relative to its installation prefix.
install: all
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/share/altann/altann $(DESTDIR)$(PREFIX)/share/altann/vendor
	install -m 755 bin/altann build/altann-core $(DESTDIR)$(PREFIX)/bin/
	install -m 644 altann/*.py $(DESTDIR)$(PREFIX)/share/altann/altann/
	install -m 644 vendor/*.pl $(DESTDIR)$(PREFIX)/share/altann/vendor/
	install -m 644 LICENSE THIRD_PARTY.md $(DESTDIR)$(PREFIX)/share/altann/

clean:
	rm -rf build

-include $(OBJECTS:.o=.d)
