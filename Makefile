# Copyright (C) 2020 University of Southern California
#                    Andrew D. Smith and Liz Ji
#
# Authors: Andrew D. Smith
#
# This file is part of EPIEVO.
#
# EPIEVO is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EPIEVO is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

all:
	@make -C src

install:
	@make -C src install

clean:
	@make -C src clean
.PHONY: clean

distclean: clean
	@rm -rf $(EPIEVO)/bin
.PHONY: distclean
