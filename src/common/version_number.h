// Copyright (C) 2012  Andrew H. Chan, Paul A. Jenkins, Yun S. Song
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// <http://www.gnu.org/licenses/>
//
// 

#ifndef LDHELMET_COMMON_VERSION_NUMBER_H_
#define LDHELMET_COMMON_VERSION_NUMBER_H_

// *VERSION_NUMBER : integer-encoding for version number (version number * 10)
// *VERSION_STRING : string of version number appropriate for display
// *_VERSION_OUTPUT: version of file output format

#define VERSION_NUMBER              19L
#define VERSION_STRING              "1.9"

#define FIND_CONFS_VERSION          10L
#define FIND_CONFS_VERSION_STRING   "1.0"

#define TABLE_GEN_VERSION            11L
#define TABLE_GEN_VERSION_STRING     "1.1"
#define TABLE_GEN_VERSION_OUTPUT     11L

#define PADE_VERSION                10L
#define PADE_VERSION_STRING         "1.0"
#define PADE_VERSION_OUTPUT         10L

#define RJMCMC_VERSION              10L
#define RJMCMC_VERSION_STRING       "1.0"
#define RJMCMC_VERSION_OUTPUT       10L

#define POST_TO_TEXT_VERSION        10L
#define POST_TO_TEXT_VERSION_STRING "1.0"

#define MAXLK_VERSION               10L
#define MAXLK_VERSION_STRING        "1.0"

#define TABLE_GEN_SALT 0xdfdf000000000000L
#define PADE_SALT      0xefef000000000000L
#define RJMCMC_SALT    0xbfbf000000000000L

#endif  // LDHELMET_COMMON_VERSION_NUMBER_H_
