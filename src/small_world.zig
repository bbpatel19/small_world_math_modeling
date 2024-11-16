const std = @import("std");

// for small world initialization
pub const n = 2000; // number of people
pub const k = 10; // number of neighbors on each side (k/2)
pub const b = 0.10; // rewire prob

var rng = std.Random.Pcg.init(1); // seed

pub const Node = u16;
pub const Edge = [2]Node; // preferrably always the right orientation so that [i,j] j-i < k
pub const Graph = [n * k]Edge;
pub const AdjMatrix = [n][n]Node;

/// the distance between two nodes on the circle, purely going right
/// assumed i,j are valid nodes
/// right_distance
pub inline fn rdist(comptime T: type, i: T, j: T) T {
    return if (j >= i) j - i else j + n - i;
}

/// starting at node i walk right j steps
/// assumed i,j are valid nodes
/// right_walk
inline fn rwalk(i: Node, j: Node) Node {
    return rdist(Node, n, i + j);
}

/// generates the small world network desired
// returns both the collection of edges in the graph and the adjacency matrix
pub fn generateSmallWorld(edges: *Graph, matrix: *AdjMatrix) void {
    {
        var in: usize = 0;
        var i: Node = 0;
        while (i < n) : (i += 1) {
            var j: Node = rwalk(i, 1);
            const end: Node = rwalk(j, k);
            while (j != end) : ({
                j = rwalk(j, 1);
            }) {
                edges[in] = .{ i, j };
                in += 1;
            }
        }
    } // this is init of graph
    for (0..n) |i| {
        for (0..i) |j| {
            const result1 = rdist(Node, @intCast(j), @intCast(i)) <= k;
            const result2 = rdist(Node, @intCast(j), @intCast(i)) >= n - k;
            matrix.*[i][j] = @intFromBool(!(result1 or result2));
            matrix.*[j][i] = @intFromBool(!(result1 or result2));
        }
        matrix.*[i][i] = 0;
    } // this is init of matrix (we keep track of this for efficiency)
    // to be noted this is actually inverted

    for (edges) |*edge| {
        if (rng.random().float(f32) < b) { // you have been chosen to be rewired
            // find new node to connect to
            const i: usize = @intCast(edge.*[0]);
            const j: usize = @intCast(edge.*[1]);
            const nn = rng.random().weightedIndex(Node, &matrix.*[i]);

            edge.* = .{ edge.*[0], @intCast(nn) };
            matrix.*[i][nn] = 0;
            matrix.*[nn][i] = 0;
            matrix.*[i][j] = 1;
            matrix.*[j][i] = 1;
        }
    }

    for (0..n) |i| { // invert the matrix so its the righ way around
        for (0..i) |j| {
            matrix.*[i][j] = 1 - matrix.*[i][j];
            matrix.*[j][i] = 1 - matrix.*[j][i];
        }
    }
}
