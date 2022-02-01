#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a neo.in file to a fortran namelist (file).

This function expects the a 'neo.in' file as input (as a string), and
will converts it to a fortran namelist 'neoin' (as string).
Comments in the code are retained.
There is also a function to read in a given file, and append the
resulting namelist to a fortran namelist file.

Note: the functions use text operations, so far no use is made of the
  f90nml class.

This file can also be used as script, in this case the append operation
is done with default file names.
"""

import sys


def all_numbers(list_input: list):
  """Check if all elements of a list can be converted to integer.

  This function expects a list and will try to convert each element in the
  list to an integer. If this is possible it will return True, and False
  otherwise.
  Note that a False is returned as soon as an element can not be
  converted, i.e. this corresponds to a short circuited evaluation.

  input:
  ------
  list_input: list, for which entries should be checked.

  output:
  -------
  logical: True if all elements can be converted to an integer, False
    otherwise.

  sideeffects:
  None.
  """
  for element in list_input:
    try:
      int(element)
    except ValueError:
      return False

  return True


def isstring(to_check):
  """Check if a given input is a string.

  Assumed to be a string, if casts to int and float fail.

  input:
  ------
  to_check: variable of which the type should be checked.

  output:
  -------
  logical: True if conversions to int and float fail, otherwise False.

  sideeffects:
  ------------
  None
  """
  try:
    int(to_check)
    return False
  except ValueError:
    try:
      float(to_check)
      return False
    except:
      return True

  return False


def convert_neo_in_to_fortran_namelist(filetext):
  """Convert given text to a fortran namelist file text.

  Comments are kept in the conversion process.

  input:
  ------
  filetext: list of strings, the file which should be converted.

  output:
  -------
  list of strings, forming a text representation of the corresponding
    namelist.
  """
  listindex=0
  outputtext = "&neoin\n"
  temp_line = ''
  arrayname = ['fluxs_arr', 'li_minima']

  for line in filetext.splitlines(keepends=True):
    # Skip empty lines, i.e. those that contain nothing, or only
    # whitespace. For safety reasons these two cases are checked
    # independently.
    if len(line) > 0:
      if (not line.isspace()):
        # The '#' is the comment character in the file.
        if line[0] == "#":
          # But fortran uses '!', so it has to be replaced.
          temp_line = None
          temp_line = '!' + line[1:]
          outputtext += temp_line
        else:
          words = line.split()
          # Check if the line contains just number, i.e. if it is an
          # array. This is the case if there are either just numbers, or
          # if the first part with a split according to '-' contains
          # just numbers. The first element for this second test has to
          # be ommited, because it might contain a minus sign.
          if all_numbers(words) or (len(line[2:].split('-')) > 1 and all_numbers(line.split('-')[0].split())):
            outputtext += " " + arrayname[listindex] + " = "
            temp = ""
            rest = ""
            for element in words:
              try:
                int(element)
              except:
                rest += ' ' + element
                continue

              temp += element + ", "

            temp += '!' + rest
            outputtext += temp + "\n"
            listindex += 1
          else:
            # Check if the value is text, i.e. a string. If so it must
            # be escaped.
            if isstring(words[0]):
              temp = " " + words[1] + " = '" + words[0] + "'   !"
            else:
              temp = " " + words[1] + " = " + words[0] + "   !"

            # Add the remainder of the line. Note that join uses the
            # string it is attached to as seperator, i.e. here it puts
            # ' ' between the words. If this would not be done we would
            # just get the text without the whitespace.
            outputtext += temp + " ".join(words[2:]) + "\n"

  outputtext += "\n/\n"

  return outputtext


def append_neo_in_to_fortran_input_file(neoinfilename: str, fortranfilename: str):
  """Convert neo.in file, and append fortran namelist file.

  input:
  ------
  neoinfilename: string, name of the neo.in file which to convert.
  fortranfilename: name of the fortran namelist file, to which to append
    the new namelist.

  output:
  -------
  None

  sideeffects:
  ------------
  Changes the file 'fortranfilename'.
  """
  f = open(neoinfilename, "r")
  filetext = f.read().strip()
  f.close()

  outputtext = convert_neo_in_to_fortran_namelist(filetext)

  # As we can not read and write at the same time(?), first read the
  # file, then reopen it to append the new namelist.
  f = open(fortranfilename, "r")
  filetext = f.read().strip()
  f.close()

  g = open(fortranfilename, "a")

  g.write(outputtext)

  g.close()


if __name__ == "__main__":
  # If this is run as main, just run the main function with the default
  # names for the files.
  append_neo_in_to_fortran_input_file('neo.in', 'neo2.in')
