const std = @import("std");
const small_world = @import("./small_world.zig");

const Node = small_world.Node;
const Graph = small_world.Graph;
const AdjMatrix = small_world.AdjMatrix;

var rng = std.Random.Pcg.init(1); // seed, what's that

// for SIR simulating
const n = small_world.n; // the total number of nodes
const initial_infected = 0.0228; // the initial infected proportion

const infection_rate = 0.0807; // the probability (weighted for each connected infected node)
const recovery_rate = 0.0415; // the probability to be recovered
// the above probability are assumed as per day
// we then approximate the prob for each time step by doing rate / num_time_steps_day
const day_infection_rate: f32 = infection_rate / @as(f32, @floatFromInt(num_time_steps_day));
const day_recovery_rate: f32 = recovery_rate / @as(f32, @floatFromInt(num_time_steps_day));

const num_time_steps_day = 100; // the number of time steps to run
const num_days = 90; // the number of days to simulate, 3 months
const iter = 10; // number of iterations to average over

const quarantine = 0.10; // prof of removing someone a distance dist-scale away
const dist_scale: f32 = 10; // distance scaling for prob of removal

const Status = enum {
    susceptible,
    infected,
    recovered,
};

fn countPopulation(nodes: [n]Status, counts_SIR: *[3]u32) void {
    counts_SIR.* = @splat(0);
    for (nodes) |node| {
        switch (node) {
            .susceptible => {
                counts_SIR.*[0] += 1;
            },
            .infected => {
                counts_SIR.*[1] += 1;
            },
            .recovered => {
                counts_SIR.*[2] += 1;
            },
        }
    }
}

fn gnuplotCounts(comptime length: usize, counts: *[length][3]f32) !void {
    const file = try std.fs.cwd().createFile("sir_counts.txt", .{ .truncate = true });
    defer file.close(); // and then close file
    var writer = std.io.bufferedWriter(file.writer());
    try std.fmt.format(writer.writer(), "#\tt\tS\tI\tR\n", .{});
    for (counts, 0..) |counts_SIR, time_days| {
        try std.fmt.format(writer.writer(), "\t{d:0>3.3}\t{d}\t{d}\t{d}\n", .{ time_days + 1, counts_SIR[0], counts_SIR[1], counts_SIR[2] });
    }
    try writer.flush();
}

inline fn min(i: usize, j: usize) usize {
    return if (i > j) j else i;
}

fn reduceConnections(matrix: *AdjMatrix, node: usize) void {
    for (matrix.*[node], 0..) |ele, i| {
        if (ele == 1) { // there is a connection
            const k = @as(f32, @floatFromInt(min(small_world.rdist(usize, node, i), small_world.rdist(usize, i, node)))) / dist_scale;
            const cutoff = 1 - std.math.pow(f32, 1 - quarantine, k);
            if (rng.random().float(f32) < cutoff) { // remove the node
                matrix.*[node][i] = 2;
                matrix.*[i][node] = 2;
            }
        }
    }
}

pub fn main() !void {
    std.debug.print("number: {d}, edges: {d}, rewire: {d}, initial infected: {d}, infection rate: {d}, recovery rate: {d}\n", .{ n, small_world.k * 2, small_world.b, initial_infected, infection_rate, recovery_rate });
    var matrix: AdjMatrix = undefined;
    var edges: Graph = undefined; // we just keep it just cause :)
    small_world.generateSmallWorld(&edges, &matrix);
    {
        var iter_num: u32 = 0;
        var tot_count_sir: [num_days][3]u32 = @splat(@splat(0));

        while (iter_num < iter) : (iter_num += 1) { // start the time step loop
            var nodes: [matrix.len]Status = @splat(.susceptible);

            for (&nodes, 0..) |*ele, i| {
                if (rng.random().float(f32) < initial_infected) {
                    ele.* = .infected;
                    if (quarantine > 0) reduceConnections(&matrix, i); // remove connections
                }
            }

            var counts_SIR: [3]u32 = @splat(0);
            countPopulation(nodes, &counts_SIR);
            std.debug.print("time step 000, S:{d}, I:{d}, R:{d}\n", .{ counts_SIR[0], counts_SIR[1], counts_SIR[2] });

            var time_days: u32 = 0;
            while (time_days < num_days) : (time_days += 1) {
                var time_step: u32 = 0;
                while (time_step < num_time_steps_day) : (time_step += 1) {
                    // SI -> II step
                    for (edges) |edge| {
                        const i: usize = @intCast(edge[0]);
                        const j: usize = @intCast(edge[1]);
                        if (matrix[i][j] != 1) continue; // will be useful when i start removing edges (not longer a node here)
                        if (nodes[i] == .infected and nodes[j] == .susceptible) {
                            if (rng.random().float(f32) < day_infection_rate) {
                                nodes[j] = .infected;
                                if (quarantine > 0) reduceConnections(&matrix, j); // remove connections
                                counts_SIR[0] -= 1;
                                counts_SIR[1] += 1;
                            }
                        }
                        if (nodes[j] == .infected and nodes[i] == .susceptible) {
                            if (rng.random().float(f32) < day_infection_rate) {
                                nodes[i] = .infected;
                                if (quarantine > 0) reduceConnections(&matrix, i); // remove connections
                                counts_SIR[0] -= 1;
                                counts_SIR[1] += 1;
                            }
                        }
                    }
                    // I -> R step
                    for (&nodes) |*node| {
                        if (node.* != .infected) continue;
                        if (rng.random().float(f32) < day_recovery_rate) {
                            node.* = .recovered;
                            counts_SIR[1] -= 1;
                            counts_SIR[2] += 1;
                        }
                    }
                }
                std.debug.print("time step {d:0>3}, S:{d}, I:{d}, R:{d}\n", .{ time_days + 1, counts_SIR[0], counts_SIR[1], counts_SIR[2] });
                const days_index: usize = @intCast(time_days);
                tot_count_sir[days_index][0] += counts_SIR[0];
                tot_count_sir[days_index][1] += counts_SIR[1];
                tot_count_sir[days_index][2] += counts_SIR[2];
            }
        }

        var avg_count_sir: [num_days][3]f32 = @splat(@splat(0));
        for (&avg_count_sir, tot_count_sir) |*row, row2| {
            for (row, row2) |*ele, ele2| {
                ele.* = @as(f32, @floatFromInt(ele2)) / iter;
            }
        }
        for (0..num_days) |i| {
            std.debug.print("time step{d}, average S: {d}, average I: {d}, average R: {d}\n", .{
                i + 1,
                avg_count_sir[i][0],
                avg_count_sir[i][1],
                avg_count_sir[i][2],
            });
        }
        try gnuplotCounts(num_days, &avg_count_sir);
    }
}
