# Basic while loop
./test 0.05 0.05 100.0 100.0 10 10 1000 50 25   # neutral case
./test 0.05 0.025 100.0 100.0 10 10 1000 50 25   # r selection case
./test 0.05 0.05 50.0 100.0 10 10 1000 50 25   # k selection case
./test 0.05 0.025 50.0 100.0 10 10 1000 50 25
#
# Scaling up K
#./test 0.15 0.15 150.0 300.0 10 10 500 100 25   # k selection case
#./test 0.15 0.15 200.0 400.0 10 10 500 100 25   # k selection case
#./test 0.15 0.15 250.0 500.0 10 10 500 100 25   # k selection case
#
# Increase r
#./test 0.3 0.1 150.0 150.0 10 10 500 100 25   # r selection case
#./test 0.4 0.1 150.0 150.0 10 10 500 100 25   # r selection case
#./test 0.5 0.1 150.0 150.0 10 10 500 100 25   # r selection case