# Chem279_Final_Project


1. compile

 g++ main.cpp periodic_solid.cpp -o eht_solver \
    -I/opt/homebrew/include \
    -L/opt/homebrew/lib \
    -larmadillo \
    -std=c++17

2. run

./eht_solver Input/graphene_clean.in
