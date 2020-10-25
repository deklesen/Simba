#[macro_use]
extern crate log;
extern crate simple_logger;
extern crate stopwatch;
extern crate rand;

use rand::{thread_rng, Rng};
use rand::distributions::{Exp, IndependentSample};

use std::cmp::Ordering;
use std::collections::BinaryHeap;
//use std::collections::HashSet;
use std::iter::FromIterator;

use std::collections::HashMap;


use std::env;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufWriter, Write};

use stopwatch::Stopwatch;

type Node = usize;

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum State {
    Infected,
    Susceptible,
    Recovered,
}
const RECOVERY_RATE: f64 = 1.0;     // Rate I+S->I+I
//const INFECTION_RATE: f64 = 1.7;    // Rate I->R  // vice versa of course
const SUSCEPTIBLE_RATE: f64 = 0.00000003;  // Rate R->S  // this is anyway hard coded to not happen
const HORIZON: f64 = 10000.0; // upper bound should stop before
const SAVEINTERVAL: usize = 10000;

#[derive(Debug, Copy, Clone)]
struct NodeInfo {
    state: State,
    recovery_time: f64, // only valid if state is infected
    susceptible_time: f64, // only valid if state is recovered or infected
    degree: usize,
}

#[derive(Debug, Copy, Clone)]
struct CountsAtTime {
    infected_count: usize,
    susceptible_count: usize,
    recovered_count: usize,
    current_time: f64,
    simulation_run: usize,
}

type Summary = Vec<CountsAtTime>;
type Node2Nodeinfo = Vec<NodeInfo>;
type GraphMap = Vec<Vec<Node>>; // Node to Neighbors

// TODO: recoery raus, time schritte anpassen, edgelist schreiben

#[derive(PartialEq, Debug, Clone)]
struct Event {
    value: f64,
    src_node: Node,
    target_node: Node,
    src_state: State,
    old_target_state: State,
    new_target_state: State,
}

type EventQueue = BinaryHeap<Event>;

impl Eq for Event {}
impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.value.partial_cmp(&self.value)
    }
}
impl Ord for Event {
    fn cmp(&self, other: &Event) -> Ordering {
        let ord = self.partial_cmp(other).unwrap();
        match ord {
            Ordering::Greater => Ordering::Less,
            Ordering::Less => Ordering::Greater,
            Ordering::Equal => ord,
        }
    }
}

fn str_2_state(state: String) -> State {
    if state == "S" {
        return State::Susceptible;
    } else if state == "I" {
        return State::Infected;
    } else if state == "R" {
        return State::Recovered;
    } else {
        panic!{"unkown state {}", state};
    }
}



fn draw_exp(rate: f64) -> f64{
    let exp = Exp::new(rate as f64);
    let v = exp.ind_sample(&mut rand::thread_rng());
    return v as f64;
}

fn sample_recovery_time() -> f64 {
    return draw_exp(RECOVERY_RATE);
}

fn sample_susceptible_time() -> f64 {
    //return draw_exp(SUSCEPTIBLE_RATE);  // change this to infinity maybe
    return 100000000.0;
}

fn get_random_neighbor(graph: &GraphMap, n: usize) -> usize {
    let neighbors = &graph[n];
    let mut rng = rand::thread_rng();
    let sample = rng.gen_range(0, neighbors.len());
    return neighbors[sample];
}

fn sample_infaction_time(degree: usize, INFECTION_RATE: f64) -> f64 {
    return draw_exp(degree as f64 * INFECTION_RATE);
}

fn infection_applicable(infection_time: f64, node: usize, node_2_nodeinfo: &Node2Nodeinfo) -> bool {
    if node_2_nodeinfo[node].state == State::Susceptible {
        return true;
    }
    if node_2_nodeinfo[node].susceptible_time < infection_time {
        return true;
    }
    return false;
}

fn create_infection_event(
    node: usize,
    node_2_nodeinfo: &Node2Nodeinfo,
    event_queue: &mut EventQueue,
    graph: &GraphMap,
    current_time: f64,
    infection_rate: f64
) {
    let mut application_time = current_time;
    loop {
        // to avoid double attack modelling
        assert!(node_2_nodeinfo[node].recovery_time >current_time);
        application_time += sample_infaction_time(node_2_nodeinfo[node].degree, infection_rate);
        if application_time>node_2_nodeinfo[node].recovery_time {
            return;
        }

        let neighbor = get_random_neighbor(graph,node);
        if infection_applicable(application_time, neighbor, node_2_nodeinfo) {
            let inf_event = Event {
                value: application_time,
                src_node: node,
                target_node: neighbor,
                src_state: State::Infected,
                old_target_state: State::Susceptible,
                new_target_state: State::Infected,
            };
            event_queue.push(inf_event);
            return;
        }
    }
}

fn setup_infection_times(node_2_nodeinfo: &mut Node2Nodeinfo, event_queue: &mut EventQueue, graph: &GraphMap, infection_rate: f64) {
    for n in 0..node_2_nodeinfo.len() {
        //let mut node_info = &mut node_2_nodeinfo[n];
        let mut waiting_time = 0.0;
        if node_2_nodeinfo[n].state == State::Infected {
            create_infection_event(n, node_2_nodeinfo, event_queue, graph, 0.0, infection_rate);
        }
    }
}

fn create_recovery_and_susceptible_event(
    node: usize,
    node_2_nodeinfo: &mut Node2Nodeinfo,
    event_queue: &mut EventQueue,
    current_time: f64,
) {
    let recovery_time = sample_recovery_time() + current_time;
    let susceptible_time = sample_susceptible_time() + recovery_time;
    let node_info = &mut node_2_nodeinfo[node];
    assert!(node_info.state == State::Infected);
    node_info.recovery_time = recovery_time;
    node_info.susceptible_time = 1000000000.0; //susceptible_time;

    let rec_event = Event {
        value: recovery_time,
        src_node: node,
        target_node: node,
        src_state: State::Infected,
        old_target_state: State::Infected,
        new_target_state: State::Recovered,
    };
    event_queue.push(rec_event);
    let sus_event = Event { // create event but do not use
        value: susceptible_time,
        src_node: node,
        target_node: node,
        src_state: State::Recovered,
        old_target_state: State::Recovered,
        new_target_state: State::Susceptible,
    };
    //event_queue.push(sus_event);
}

fn setup_recovery_times(node_2_nodeinfo: &mut Node2Nodeinfo, event_queue: &mut EventQueue) {
    for n in 0..node_2_nodeinfo.len() {
        if node_2_nodeinfo[n].state == State::Infected {
            create_recovery_and_susceptible_event(n, node_2_nodeinfo, event_queue, 0.0);
        }
    }
}

fn setup_graph(
    graphpath: String,
    graph: &mut GraphMap,
    node_infos: &mut Node2Nodeinfo,
    current_counts: &mut CountsAtTime,
    initially_infected: &mut Vec<usize>,
) {
    let mut file = match File::open(&graphpath) {
        Err(_why) => panic!("couldn't find graphfile"),
        Ok(file) => file,
    };
    let mut s = String::new();
    match file.read_to_string(&mut s) {
        Err(_why) => panic!("couldn't read graphfile"),
        Ok(_) => (),
    }
    let lines = s.split("\n");
    let mut label: State;
    let mut counter: usize = 0;
    let mut degree: usize;

    for l in lines {
        if l.len() < 3 {
            continue;
        }
        let line_info: Vec<&str> = l.split(";").collect();

        if line_info[0].to_string().parse::<usize>().unwrap() != counter {
            println!("Wrong order of nodes in input graph");
        }
        counter += 1;
        degree = 0;

        //let v: &str = line_info[0];
        //v = line_info[0].to_string().parse().unwrap();   TODO why is this never used
        label = str_2_state(line_info[1].to_string());
        if label == State::Infected {
            current_counts.infected_count += 1;
            //println!("Initially infected {:?}", graph.len());
            initially_infected.push(graph.len());
        } else if label == State::Susceptible {
            current_counts.susceptible_count += 1;
        } else {
            current_counts.recovered_count += 1;
        }
        if line_info[2].len() > 0 {
            // TODO check if only one neighbor
            let neighbors_str: Vec<&str> = line_info[2].split(",").collect();
            degree = neighbors_str.len();
            let neighbors: Vec<Node> = neighbors_str
                .iter()
                .map(|v| v.to_string().parse::<Node>().unwrap())
                .collect();
            graph.push(neighbors);
            node_infos.push(NodeInfo {
                state: label,
                recovery_time: 0.0,
                susceptible_time: 0.0,
                degree: degree,
            });
        } else {
            println!("Node without neighbour occured{:?}", l);
            //panic!("use kmax larger than 0");
            let neighbors: Vec<Node> = [].to_vec();
            graph.push(neighbors);
            node_infos.push(NodeInfo {
                state: label,
                recovery_time: 0.0,
                susceptible_time: 0.0,
                degree: degree,
            });
        }
    }
}

fn read_arguments() -> (String, String, String, f64, usize, String) {
    //(graphpath, outpath_traj, outpath_TG, outpath_score, infection_rate)
    let args: Vec<String> = env::args().collect();
    return (args[1].clone(), args[2].clone(), args[3].clone(), args[4].clone().to_string().parse::<f64>().unwrap(), args[5].clone().to_string().parse::<usize>().unwrap(),args[6].clone());
}



fn apply_recovery(current_event: &Event, node_2_nodeinfo: &mut Node2Nodeinfo) {
    let node = current_event.src_node;
    info!("recover! {}", node);
    assert!(node_2_nodeinfo[node].state == State::Infected);
    node_2_nodeinfo[node].state = State::Recovered;
}

fn apply_susceptible(current_event: &Event, node_2_nodeinfo: &mut Node2Nodeinfo) {
    let node = current_event.src_node;
    info!("susceptible! {}", node);
    assert!(node_2_nodeinfo[node].state == State::Recovered);
    node_2_nodeinfo[node].state = State::Susceptible;
}

fn apply_infection(current_event: &Event, node_2_nodeinfo: &mut Node2Nodeinfo) -> bool {
    let src_node = current_event.src_node;
    let target_node = current_event.target_node;
    assert!(src_node != target_node);
    assert!(current_event.src_state == State::Infected);
    assert!(node_2_nodeinfo[src_node].state == State::Infected);
    assert!(current_event.old_target_state == State::Susceptible);
    if node_2_nodeinfo[src_node].state == current_event.src_state
        && node_2_nodeinfo[target_node].state == current_event.old_target_state
    {
        info!("infect! {}", target_node);
        node_2_nodeinfo[target_node].state = State::Infected;
        return true;
    }
    return false;
}

fn apply_event(current_event: &Event, node_2_nodeinfo: &mut Node2Nodeinfo) -> bool {
    // if recovery event
    if current_event.new_target_state == State::Recovered {
        apply_recovery(&current_event, node_2_nodeinfo);
        return true; // is always successful
    } else if current_event.new_target_state == State::Infected {
        // is infection event
        let was_successful = apply_infection(&current_event, node_2_nodeinfo);
        return was_successful;
    } else {
        // is infection event
        assert!(current_event.new_target_state == State::Susceptible);
        apply_susceptible(&current_event, node_2_nodeinfo);
        return true;
    }
}

fn perform_step(
    graph: &mut GraphMap,
    node_2_nodeinfo: &mut Node2Nodeinfo,
    current_counts: &mut CountsAtTime,
    event_queue: &mut EventQueue,
    current_time: f64,
    tg: &mut HashMap<(usize,usize), f64>,
    infection_rate: f64
) -> (f64, bool) {
    // only updates counts not time in current_counts

    if event_queue.len() == 0 {
        info!("No events left. Simualtion over.");
        return (HORIZON + 0.0000001, false);
    }
    let current_event = event_queue.pop().unwrap();
    let current_time = current_event.value;
    info!("{:?}", current_event);

    let was_successful = apply_event(&current_event, node_2_nodeinfo);
    info!("was successful {}", was_successful);

    //
    // TODO add check
    //
    if was_successful && current_event.new_target_state == State::Infected && current_event.target_node != current_event.src_node  {
        let mut count: f64 = 0.0;
        if tg.contains_key(&(current_event.target_node, current_event.src_node)){
            count = tg[&(current_event.target_node, current_event.src_node)];
        }
        tg.insert((current_event.target_node, current_event.src_node),count+1.0);

        count = 0.0;
        if tg.contains_key(&(current_event.target_node, current_event.target_node)){
            count = tg[&(current_event.target_node, current_event.target_node)];
        }
        tg.insert((current_event.target_node, current_event.target_node),count+1.0);

        //if !tg.contains_key(&(current_event.src_node, current_event.src_node)){
        //    tg.insert((current_event.src_node, current_event.src_node),0.0);
        //}
    }


    //current_event.is_infection
    if current_event.new_target_state == State::Infected {
        if was_successful {
            // order is important, recovery before infection event
            create_recovery_and_susceptible_event(
                current_event.target_node,
                node_2_nodeinfo,
                event_queue,
                current_time,
            );
            create_infection_event(
                current_event.target_node,
                node_2_nodeinfo,
                event_queue,
                graph,
                current_time,
                infection_rate
            );
        }
        create_infection_event(
            current_event.src_node,
            node_2_nodeinfo,
            event_queue,
            graph,
            current_time,
            infection_rate
        );
    }

    // successful infection event
    if was_successful && (current_event.new_target_state == State::Infected)  {
        current_counts.susceptible_count -= 1;
        current_counts.infected_count += 1;
    } else if was_successful && (current_event.new_target_state == State::Susceptible) {
        current_counts.susceptible_count += 1;
        current_counts.recovered_count -= 1;
    } else if was_successful && (current_event.new_target_state == State::Recovered) {
        current_counts.infected_count -= 1;
        current_counts.recovered_count += 1;
    }

    return (current_time, was_successful);
}

fn save_system_state(summary: &mut Summary, current_time: f64, current_counts: CountsAtTime) {
    summary.push(current_counts);
    return;
}

fn write_output(summary: Summary, step_count: usize, rejected_steps: usize, runtime: i64, outpath: String) {
    let outpath_runtime = outpath.replace(".txt", "_runtime.txt");

    let output_rows = 1000;
    let result_len = summary.len() as i32;
    let subsamp_index = result_len/output_rows as i32;
    let mut counter = 0;

    let mut f = BufWriter::new(fs::File::create(outpath_runtime).unwrap());
    write!(f, "runtime(ms),steps,rejected_steps\n{:?},{},{}", runtime, step_count,rejected_steps);

    let mut f = BufWriter::new(fs::File::create(outpath).unwrap());
    write!(f, "state,fraction,time\n");
    for node_count in summary {
        counter += 1;
        if result_len > output_rows*2 && counter > 100 && counter < result_len-100 && counter % subsamp_index != 0  {
            continue;
        }

        let number_of_nodes = (node_count.susceptible_count + node_count.infected_count + node_count.recovered_count) as f32;
        let s_frac = (node_count.susceptible_count as f32) / number_of_nodes;
        let i_frac = (node_count.infected_count as f32) / number_of_nodes;
        let r_frac = (node_count.recovered_count as f32) / number_of_nodes;
        let time = node_count.current_time;
        write!(f, "S,{},{}\n", s_frac, time);
        write!(f, "R,{},{}\n", r_frac, time);
        write!(f, "I,{},{}\n", i_frac, time);
    }

    return;
}

fn print_event_queue(mut event_queue: EventQueue) {
    while !event_queue.is_empty(){
        println!("{:?}", event_queue.pop());
    }
}


fn main_single(mut graph: GraphMap, mut current_counts: CountsAtTime, mut node_2_nodeinfo: Node2Nodeinfo, tg: &mut HashMap<(usize,usize), f64>, mut summary: &mut  Summary, step_i: usize, infection_rate: f64) -> f64 {

    print!(".");
    let mut event_queue: EventQueue = BinaryHeap::new();

    setup_recovery_times(&mut node_2_nodeinfo, &mut event_queue);
    setup_infection_times(&mut node_2_nodeinfo, &mut event_queue, &graph, infection_rate);
    info!("current counts {}", current_counts.susceptible_count);

    //print_event_queue(event_queue.clone());

    let mut current_step: usize = 0;
    let mut current_time: f64 = 0.0;
    let mut real_steps: usize = 0;
    let mut rejected_steps: usize = 0;


    while current_time < HORIZON {
        info!("Current counts: {:?}, time: {}", current_counts, current_time);

        if current_step < 100 || current_step % SAVEINTERVAL == 0 {
            current_counts.simulation_run = step_i;
            save_system_state(summary, current_time, current_counts.clone());
        }

        let (event_time, was_successful) = perform_step(
            &mut graph,
            &mut node_2_nodeinfo,
            &mut current_counts,
            &mut event_queue,
            current_time,
            tg,
            infection_rate
        );
        current_time = event_time;
        current_step += 1;
        if was_successful {
            real_steps += 1;
        } else {
            rejected_steps += 1;
        }

        current_counts.current_time = current_time;

        if current_step % 10000 == 0 {
            //print!(".");
            if current_step % 1000000 == 0 {
                println!("\ntime: {}", current_time);
            }
        }
    }
     return (current_counts.susceptible_count as f64)/((current_counts.susceptible_count as f64)+(current_counts.recovered_count as f64)+(current_counts.infected_count as f64));
}




fn write_summary(summary: Summary, outpath: &String) {

    let output_rows = 1000;
    let result_len = summary.len() as i32;
    let subsamp_index = result_len/output_rows as i32;
    let mut counter = 0;

    let mut f = BufWriter::new(fs::File::create(outpath).unwrap());
    write!(f, "State,Fraction,Time,Simulation_run\n");
    for node_count in summary {
        counter += 1;
        if result_len > output_rows*2 && counter > 100 && counter < result_len-100 && counter % subsamp_index != 0  {
            continue;
        }

        let number_of_nodes = (node_count.susceptible_count + node_count.infected_count + node_count.recovered_count) as f32;
        let s_frac = (node_count.susceptible_count as f32) / number_of_nodes;
        let i_frac = (node_count.infected_count as f32) / number_of_nodes;
        let r_frac = (node_count.recovered_count as f32) / number_of_nodes;
        let time = node_count.current_time;
        write!(f, "S,{},{},{}\n", s_frac, time, node_count.simulation_run);
        write!(f, "R,{},{},{}\n", r_frac, time, node_count.simulation_run);
        write!(f, "I,{},{},{}\n", i_frac, time, node_count.simulation_run);
    }

    return;
}


// in rust_reject  cargo build --release
// call from root with ./rust_reject/target/release/rust_reject example_networks/chain_graph.txt out_trajectory.txt out_TG 2.2 33        //
fn main() {
    simple_logger::init_with_level(log::Level::Warn).unwrap();
    let (graphpath, outpath, outpath_TG, infection_rate, num_run, write_output_string) = read_arguments();
    let write_output: bool = (write_output_string=="Yes");
    println!("Input: {:?}", (&graphpath, &outpath, &outpath_TG, infection_rate, num_run, write_output, write_output_string));
    let stopwatch = Stopwatch::start_new();

    let mut current_counts: CountsAtTime = CountsAtTime {
    infected_count: 0,
    recovered_count: 0,
    susceptible_count: 0,
    current_time: 0.0,
    simulation_run: 0,
    };
    let mut node_2_nodeinfo: Node2Nodeinfo = Vec::new();

    let mut summary: Summary = Vec::with_capacity(1000000);

    let mut initially_infected: Vec<usize> = Vec::new();
    let mut graph: GraphMap = Vec::new();
    setup_graph(
        graphpath,
        &mut graph,
        &mut node_2_nodeinfo,
        &mut current_counts,
        &mut initially_infected,
    );

    let mut transmission_graph:HashMap<(usize,usize), f64> = HashMap::new();
    for i in 0..graph.len() {
        let v = &graph[i];
        for j in 0..v.len() {
            if !initially_infected.contains(&i){
            transmission_graph.insert((i,v[j]), 0.001);
            transmission_graph.insert((i,i), 0.001);
            }
        }
    }

    let mut transmission_graph_normalized:HashMap<(usize,usize), f64> = HashMap::new();
    //transmission_graph.insert((0,0), 0.0);
    let mut final_unaffected: Vec<f64> = Vec::new();
    let mut score: f64 = 0.0;

    for step_i in 0..num_run {
        score = main_single(graph.clone(), current_counts.clone(), node_2_nodeinfo.clone(), &mut transmission_graph, &mut summary, step_i, infection_rate);
        final_unaffected.push(score);
    }

    let elapsed_time = stopwatch.elapsed_ms();
    println!("Elapsed Time: {:?} ms", &elapsed_time);

    //for (edge, counts) in &transmission_graph {
    //    println!("{:?}: \"{:?}\"", edge, counts);
    //}



    if write_output{
            
        let mut f = BufWriter::new(fs::File::create(&outpath_TG).unwrap());
        let mut normalize = 0.0;
        for (edge, counts) in &transmission_graph {
            if transmission_graph.contains_key(&(edge.0, edge.0)) && edge.0 != edge.1{
                normalize = transmission_graph[&(edge.0, edge.0)];
                if normalize > 0.0 {
                    transmission_graph_normalized.insert((edge.0, edge.1), counts/normalize);
                    write!(f, "{:?} {:?} {:?}\n", edge.0, edge.1, counts/normalize);
                }
            }
        }

        let mut inf_total = 0.0;
        for node in 0..graph.len() {
            if transmission_graph.contains_key(&(node, node)) {
                inf_total = inf_total +  transmission_graph[&(node, node)]
            }
        }
        let dummy_node: usize = graph.len();

        for node in 0..graph.len() {
            if transmission_graph.contains_key(&(node, node)) {
                let prob = transmission_graph[&(node,node)] / inf_total;
                transmission_graph_normalized.insert((dummy_node, node), prob);
                write!(f, "{:?} {:?} {:?}\n", dummy_node, node, prob);
            }
        }

        for node in &initially_infected {
            transmission_graph_normalized.insert((*node, dummy_node), 1.0);
            write!(f, "{:?} {:?} {:?}\n", node, dummy_node, 1.0);
        }

        println!("len summary: {:?}", summary.len());
        write_summary(summary, &outpath);

    }

    // write score:
    let mut f = BufWriter::new(fs::File::create(format!("{}.score", &outpath)).unwrap());
    for score in &final_unaffected {
        write!(f, "{:?}\n", score);
    }

    if write_output{
        // write intensity (numer of times nodes became infected):
        // is zero for initially infected nodes
        let mut f = BufWriter::new(fs::File::create(format!("{}.intensity", &outpath_TG)).unwrap());
        for node in 0..graph.len() {
            if transmission_graph.contains_key(&(node, node)) {
                write!(f, "{:?}\n", transmission_graph[&(node, node)]);
            } else {
                write!(f, "0\n");
            }
        }
    }


    //for (edge, counts) in &transmission_graph_normalized {
    //    println!("{:?}: \"{:?}\"", edge, counts);
    //}
    if write_output{
        // solve:  TG
        let number_of_nodes = graph.len() + 1; // +1 is dummy
        let mut p_dist : Vec<f64> = vec![1.0/(number_of_nodes as f64); number_of_nodes];
        for step_i in 0..1000 {
            let mut p_dist_new : Vec<f64> = vec![0.0; number_of_nodes];
            for (edge, weight) in &transmission_graph_normalized {
                if edge.0 != edge.1 {
                    p_dist_new[edge.1] += p_dist[edge.0] * weight;
                }
            }
            p_dist = p_dist_new.clone();

            // normalize against numerical problems
            let normalize: f64 = p_dist.iter().sum();
            for i in 0..p_dist.len(){
                p_dist[i] = p_dist[i]/normalize;
            }
        }

        let lel = p_dist.len()-1;
        p_dist[lel] = 0.0;
        for node in initially_infected {
            p_dist[node] = 0.0;
        }
        let normalize: f64 = p_dist.iter().sum();
        for i in 0..p_dist.len(){
            p_dist[i] = p_dist[i]/normalize;
        }

        let mut f = BufWriter::new(fs::File::create(format!("{}.solution", &outpath_TG)).unwrap());
        for i in 0..p_dist.len() {
            write!(f, "{:?}\n", p_dist[i]);
        }

        let mut f = BufWriter::new(fs::File::create(format!("{}.elapsedTime", &outpath_TG)).unwrap());
            write!(f, "Elapsed Time: {:?} ms", &elapsed_time);

    }
    println!("rust done\n");

    //println!("{:?}\n", &p_dist);


    //info!("Number of steps: {}", current_step);
    //save_system_state(&mut summary, current_time, current_counts.clone());
    //write_output(summary, real_steps, rejected_steps, elapsed_time, outpath);

}

