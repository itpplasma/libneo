%add path of locale installation of InputFile class
addpath('~/GITHUB/libneo/matlab/InputFile');

%initialize class
gorilla = InputFile('./gorilla.inp');

%read file
gorilla.read();

%make some changes
gorilla.GORILLANML.n1 = 200;
gorilla.GORILLANML.boole_guess = false;
gorilla.GORILLANML.coord_system = 1;

%write to output
gorilla.write('./gorilla.out');