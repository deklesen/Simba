# Rejection-Based Epidemic Simulation for Complex Networks
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/gerritgr/Rejection-Based-Epidemic-Simulation.svg?token=qQ7vTmAySdBppYxywojC&branch=master)](https://travis-ci.com/gerritgr/Rejection-Based-Epidemic-Simulation)
## Overview
Implementation of Monte-Carlo (Gillespie) simulation for epidemic type processes on complex networks.

## Installation
Install Rust with:
```sh
curl https://sh.rustup.rs -sSf | sh
```
Compile the Rust code (in the rust_reject folder) with:
```sh
cd rust_reject && cargo build --release  
```

## Usage

Start the SIS simulation with
```sh
./rust_reject/target/release/rust_reject example_networks/gamma_2.0_nodes_1000.txt out_trajectory.txt
```

where the first argument is the input network and the second argument is the output-filepath.

### Network File 
The network file contains containing a labeled graph specifying the initial state, each line having the form `<Nodeid>;<Label>;<Neighbor1>,<Neighbor2>,...`
```sh
0;I;31,29,94,13,83
1;S;66,15,73
2;S;29,61,26,80,16,83,30,62,3,93,27,87,68,18,79,6
3;I;83,2,29,4,28,61,46,21,9,49,41,68,16,74
4;S;82,28,12,83,3,62,66,68
...
```
Nodes start with id 0 and are sorted. 
Isolates (nodes withouth neighbors) are not supported (yet). 
There should be at least one node for each possible label. 

### Rates
Currently, it is only possible to specify the rate parameters directly in the `main.rs` source file:
```sh
const RECOVERY_RATE: f64 = 1.0;
const INFECTION_RATE: f64 = 0.6;
const HORIZON: f64 = 10.0;
const SAVEINTERVAL: usize = 1000;
```
SAVEINTERVAL specifies the time resolution of the output. 

### SIR and Competing Pathogens Model
To use one of the other models simply rename them to `main.rs`.
