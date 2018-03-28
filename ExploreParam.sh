# Basic while loop
./test 0.1 0.1 150.0 150.0 10 10 500 100 25   # neutral case
./test 0.2 0.1 150.0 150.0 10 10 500 100 25   # r selection case
./test 0.15 0.15 100.0 200.0 10 10 500 100 25   # k selection case

# Scaling up K
./test 0.15 0.15 150.0 300.0 10 10 500 100 25   # k selection case
./test 0.15 0.15 200.0 400.0 10 10 500 100 25   # k selection case
./test 0.15 0.15 250.0 500.0 10 10 500 100 25   # k selection case

# Increase r
./test 0.3 0.1 150.0 150.0 10 10 500 100 25   # r selection case
./test 0.4 0.1 150.0 150.0 10 10 500 100 25   # r selection case
./test 0.5 0.1 150.0 150.0 10 10 500 100 25   # r selection case