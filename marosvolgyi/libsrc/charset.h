/*
characterset for libdraw.
    Copyright (C) 2010 Marcell Marosvolgyi, a.k.a. cello

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
short A[64*67] = {
  //0 0
  0,1,1,1,1,1,0,0,
  1,0,0,0,0,1,1,0,
  1,0,0,0,1,0,1,0,
  1,0,0,1,0,0,1,0,
  1,0,1,0,0,0,1,0,
  1,1,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //1 1
  0,0,1,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //2 2
  0,0,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,1,1,1,1,0,0,
  0,1,0,0,1,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0,
  //3 3
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //4 4
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //5 5
  1,1,1,1,1,1,1,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //6 6
  0,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //7 7
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //8 8
  0,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //9 9
  0,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,1,1,1,1,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  // ! 10
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // @ 11
  0,0,1,1,1,0,0,0,
  0,1,0,0,0,1,0,0,
  1,0,0,1,0,0,1,0,
  1,0,1,0,1,0,1,0,
  1,0,1,0,1,0,1,0,
  1,0,0,1,0,1,0,0,
  0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // # 12
  0,0,1,0,1,0,0,0,
  0,0,1,0,1,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,1,0,1,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,1,0,1,0,0,0,
  0,0,1,0,1,0,0,0,
  0,0,0,0,0,0,0,0,
  // $ 13
  0,0,0,1,0,0,0,0,
  0,0,1,1,1,0,0,0,
  0,1,0,1,0,0,0,0,
  0,0,1,1,0,0,0,0,
  0,0,0,1,1,0,0,0,
  0,0,1,1,0,1,0,0,
  0,0,0,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  // % 14
  0,1,0,0,0,0,1,0,
  1,0,1,0,0,1,0,0,
  0,1,0,0,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,1,0,0,0,0,1,0,
  1,0,0,0,0,1,0,1,
  0,0,0,0,0,0,1,0,
  // ^15
  0,0,0,1,0,0,0,0,
  0,0,1,0,1,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // &16
  0,0,1,1,0,0,0,0,
  0,1,0,0,1,0,0,0,
  0,0,1,0,0,0,0,0,
  0,1,0,1,0,0,1,0,
  1,0,0,0,1,1,0,0,
  1,0,0,0,1,1,0,0,
  0,1,1,1,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //* 17
  1,0,0,1,0,0,1,0,
  0,1,0,1,0,1,0,0,
  0,0,1,1,1,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,1,1,1,0,0,0,
  0,1,0,1,0,1,0,0,
  1,0,0,1,0,0,1,0,
  0,0,0,0,0,0,0,0,
  // ( 18
  0,0,1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // ) 19
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //+ 20
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //- 21
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // / 22
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // _ 23
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0,
  // < 24
  0,0,0,1,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // > 25
  0,0,0,1,1,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // ! 26
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // ? 27
  0,0,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,1,1,1,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 28
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 29
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 30
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 27
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 31
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 32
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 33
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 34
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 35
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 36
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 37
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 38
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  // .. 39
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //A 40
  0,0,0,1,0,0,0,0,
  0,0,1,0,1,0,0,0,
  0,1,0,0,0,1,0,0,
  1,1,1,1,1,1,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //B 41
  1,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //C 42
  0,0,1,1,1,1,0,0,
  0,1,0,0,0,0,1,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,1,0,0,0,0,1,0,
  0,0,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //D 30
  1,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //E 31
  1,1,1,1,1,1,1,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0,
  //F 32
  1,1,1,1,1,1,1,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,1,1,1,1,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //G 33
  0,0,1,1,1,1,0,0,
  0,1,0,0,0,0,1,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,1,1,1,0,
  0,1,0,0,0,0,1,0,
  0,0,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //H 34
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,1,1,1,1,1,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //I 35
  0,0,1,1,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //J 36
  0,0,1,1,1,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,0,1,0,0,0,
  0,1,0,0,1,0,0,0,
  0,0,1,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //K
  1,0,0,0,0,1,1,0,
  1,0,0,0,1,0,0,0,
  1,0,0,1,0,0,0,0,
  1,1,1,1,0,0,0,0,
  1,0,0,1,0,0,0,0,
  1,0,0,0,1,0,0,0,
  1,0,0,0,0,1,1,0,
  0,0,0,0,0,0,0,0,
  //L
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,1,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0,
  //M
  1,1,0,0,0,0,1,0,
  1,0,1,0,0,1,1,0,
  1,0,0,1,1,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //N
  1,1,0,0,0,0,1,0,
  1,0,1,0,0,0,1,0,
  1,0,0,1,0,0,1,0,
  1,0,0,0,1,0,1,0,
  1,0,0,0,0,1,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //O
  0,0,1,1,1,0,0,0,
  0,1,0,0,0,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,0,0,0,1,0,0,
  0,0,1,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //P
  1,1,1,1,1,0,0,0,
  1,0,0,0,0,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,1,0,0,
  1,1,1,1,1,0,0,0,
  1,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //Q
  0,0,1,1,1,0,0,0,
  0,1,0,0,0,0,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,1,0,0,1,0,
  0,1,0,0,1,1,0,0,
  0,0,1,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //R
  1,1,1,1,1,1,0,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,1,0,0,
  1,1,1,1,1,0,0,0,
  1,0,0,1,0,0,0,0,
  1,0,0,0,1,1,1,0,
  0,0,0,0,0,0,0,0,
  //S
  0,1,1,1,1,1,1,0,
  1,0,0,0,0,0,0,0,
  1,0,0,0,0,0,0,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  1,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //T
  1,1,1,1,1,1,1,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,0,
  //U
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,1,1,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //V
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,0,0,0,0,1,0,
  0,0,1,0,0,1,0,0,
  0,0,0,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //W
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,1,0,0,1,0,
  1,0,0,1,0,0,1,0,
  0,1,1,0,1,1,0,0,
  0,0,0,0,0,0,0,0,
  //X
  1,0,0,0,0,0,1,0,
  0,1,0,0,0,1,0,0,
  0,0,1,0,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,0,1,0,0,0,
  0,1,0,0,0,1,0,0,
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  //Y
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  1,0,0,0,0,0,1,0,
  0,1,0,0,0,0,1,0,
  0,0,1,1,1,1,1,0,
  0,0,0,0,0,1,0,0,
  1,1,1,1,1,0,0,0,
  0,0,0,0,0,0,0,0,
  //Z
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,1,0,0,
  0,0,0,0,1,0,0,0,
  0,0,0,1,0,0,0,0,
  0,0,1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,
  1,1,1,1,1,1,1,0,
  0,0,0,0,0,0,0,0
};