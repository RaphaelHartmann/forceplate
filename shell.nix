{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    bashInteractive
    rEnv
    rstudioEnv
    git
    rPackages.data_table
    rPackages.signal
    # ... any other dependencies
  ];

  # You can set environment variables as needed
  # environment variables that should be set when the shell is entered
  shellHook = ''
    echo "entering rEnv rstudioEnv and gsl evnironment"
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
    #  export SOME_VAR=some_value
  '';
}
