/*################################################################################
  ##
  ##   Copyright (C) 2017-2018 Keith O'Hara
  ##
  ##   This file is part of the Shortest Path C++ library (SPLib).
  ##
  ##   SPLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   SPLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * simple tictoc functionality
 */

using clocktime_t = std::chrono::time_point<std::chrono::system_clock>;
using comptime_t  = std::chrono::duration<double>;

inline
clocktime_t
tic()
{
    return std::chrono::system_clock::now();
}

inline
void
tictoc(clocktime_t time_inp)
{
    clocktime_t time_now = std::chrono::system_clock::now();

    comptime_t run_time = time_now - time_inp;

    //

    std::time_t end_time = std::chrono::system_clock::to_time_t(time_now);
        
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "total runtime: " << run_time.count() << "s\n";
}
