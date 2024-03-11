use ds::Matrix;
use std::fmt;
use std::fs;

#[derive(Debug)]
struct PData {
    x: f64,
    y: f64,
}

impl PData {
    fn new(x: f64, y: f64) -> PData {
        PData { x, y }
    }
}

#[derive(Debug)]
struct MData {
    x: usize,
    y: usize,
    e: f64,
    i: f64,
    a: f64,
}

impl MData {
    fn new(x: usize, y: usize, e: f64, i: f64, a: f64) -> MData {
        MData { x, y, e, i, a }
    }
}

#[derive(Debug)]
struct NodalLoad {
    x: f64,
    y: f64,
    m: f64,
    index: usize,
}

impl NodalLoad {
    fn new(x: f64, y: f64, m: f64, index: usize) -> NodalLoad {
        NodalLoad { x, y, m, index }
    }
}

#[derive(Debug)]
struct PointLoad {
    val: f64,
    start: usize,
    end: usize,
    offset: f64,
}

impl PointLoad {
    fn new(val: f64, start: usize, end: usize, offset: f64) -> PointLoad {
        PointLoad {
            val,
            start,
            end,
            offset,
        }
    }
}

#[derive(Debug)]
struct UDL {
    val: f64,
    start: usize,
    end: usize,
}

impl UDL {
    fn new(val: f64, start: usize, end: usize) -> UDL {
        UDL { val, start, end }
    }
}

#[derive(Debug)]
struct UVL {
    start_val: f64,
    end_val: f64,
    start: usize,
    end: usize,
}

impl UVL {
    fn new(start_val: f64, end_val: f64, start: usize, end: usize) -> UVL {
        UVL {
            start_val,
            end_val,
            start,
            end,
        }
    }
}

#[derive(Debug)]
struct SData {
    typ: String,
    index: usize,
    y_orient: bool,
}

impl SData {
    fn new(typ:String, index: usize, orientation: bool) -> SData {
        SData { typ, index, y_orient: orientation}
    }
}


#[derive(Debug, Copy, Clone, PartialEq)]
struct Point {
    x: f64,
    y: f64,
    id: usize,
}

impl Point {
    fn new(x: f64, y: f64) -> Point {
        Point { x, y, id: 0 }
    }

    fn new_from_int(x: usize, y: usize) -> Point {
        let xf64 = x as f64;
        let yf64 = y as f64;
        Point {
            x: xf64,
            y: yf64,
            id: 0,
        }
    }

    fn distance(&self, other: Point) -> f64 {
        ((self.y - other.y).powi(2) + (self.x - other.x).powi(2)).sqrt()
    }

    fn slope(&self, other: Point) -> f64 {
        (other.y - self.y) / (other.x - self.x)
    }
}

impl std::fmt::Display for Point {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {}) -> id={}", self.x, self.y, self.id)
    }
}

enum Support {
    Roller {
        id: usize,
        position: Point,
        y_orientation: bool,
    },
    Hinged {
        id: usize,
        position: Point,
    },
    Fixed {
        id: usize,
        position: Point,
    },
}

impl Support {
    fn new_roller(id: usize, position: Point, y_orientation: bool) -> Support {
        Support::Roller {
            id,
            position,
            y_orientation,
        }
    }

    fn new_hinged(id: usize, position: Point) -> Support {
        Support::Hinged { id, position }
    }

    fn new_fixed(id: usize, position: Point) -> Support {
        Support::Fixed { id, position }
    }

    fn get_fixed_dofs(&self) -> Vec<usize> {
        let mut fd: Vec<usize> = vec![];
        match self {
            Support::Roller {
                id, y_orientation, ..
            } => {
                if *y_orientation {
                    fd.extend(vec![3 * id + 1]);
                } else {
                    fd.extend(vec![3 * id]);
                }
            }
            Support::Hinged { id, .. } => {
                fd.extend(vec![3 * id, 3 * id + 1]);
            }
            Support::Fixed { id, .. } => {
                fd.extend(vec![3 * id, 3 * id + 1, 3 * id + 2]);
            }
        }
        fd
    }
}



#[derive(Clone)]
struct Member {
    node1: usize,
    node2: usize,
    point1: Point,
    point2: Point,
    e: f64,
    i: f64,
    a: f64,
}

impl Member {
    fn new(n1: usize, n2: usize, p1: Point, p2: Point) -> Member {
        Member {
            node1: n1,
            node2: n2,
            point1: p1,
            point2: p2,
            e: 1.0,
            i: 1.0,
            a: 1.0,
        }
    }

    fn new_withparams(n1: usize, n2: usize, p1: Point, p2: Point, e: f64, i: f64, a: f64) -> Member {
        Member{
            node1: n1,
            node2: n2,
            point1: p1,
            point2: p2,
            e,i,a,
        }
    }


    fn contains(&self, node: usize) -> bool {
        self.node1 == node || self.node2 == node
    }

    fn length(&self) -> f64 {
        self.point2.distance(self.point1)
    }

    fn slope(&self) -> f64 {
        self.point1.slope(self.point2)
    }

    fn stiffness(&self) -> Matrix {
        let t1 = self.e * self.a / self.length();
        let t2 = 12.0 * self.e * self.i / f64::powi(self.length(), 3);
        let t3 = 6.0 * self.e * self.i / f64::powi(self.length(), 2);
        let t4 = 4.0 * self.e * self.i / self.length();
        let t5 = 2.0 * self.e * self.i / self.length();

        Matrix::new(
            vec![
                t1, 0.0, 0.0, -t1, 0.0, 0.0, 0.0, t2, t3, 0.0, -t2, t3, 0.0, t3, t4, 0.0, -t3, t5,
                -t1, 0.0, 0.0, t1, 0.0, 0.0, 0.0, -t2, -t3, 0.0, t2, -t3, 0.0, t3, t5, 0.0, -t3,
                t4,
            ],
        )
    }

    fn transformation_matrix(&self) -> Matrix {
        let l = self.slope().atan().cos();
        let m = self.slope().atan().sin();
        Matrix::new(
            vec![
                l, m, 0.0, 0.0, 0.0, 0.0, 
                -m, l, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, l, m, 0.0,
                0.0, 0.0, 0.0, -m, l, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        )
    }

    fn stiffness_matrix(&self) -> Matrix {
        let t = self.transformation_matrix();
        t.transpose().multipy(self.stiffness()).multipy(t)

    }


}

impl std::fmt::Display for Member {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {})", self.node1, self.node2)
    }
}

struct Frame {
    nodes: Option<Vec<Point>>,
    members: Option<Vec<Member>>,
    nodes_count: usize,
    members_count: usize,
    fixed_dofs: Option<Vec<usize>>,
    load_vector: Option<Vec<f64>>,
    solution_displacement: Option<Vec<f64>>,
    solution_load : Option<Vec<f64>>,
}

impl fmt::Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Nodes:\n--------------\n")?;
        for (i, s) in self.get_nodes().iter().enumerate() {
            if i > 0 {
                write!(f, ",\n")?;
            }
            write!(f, "{}", s)?;
        }
        write!(f, "\n\nMembers: [")?;
        for (i, s) in self.get_members().iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{}", s)?;
        }
        write!(f, "],\tCount={}", self.nodes_count)
    }
}

impl Frame {
    fn new() -> Frame {
        Frame {
            nodes: None,
            members: None,
            nodes_count: 0,
            members_count: 0,
            fixed_dofs: None,
            load_vector: None,
            solution_displacement: None,
            solution_load: None,
        }
    }

    fn get_nodes(&self) -> Vec<Point> {
        if let Some(nodes) = &self.nodes {
            return nodes.to_vec();
        }
        panic!("nodes not defined")
    }

    fn get_members(&self) -> Vec<Member> {
        if let Some(members) = &self.members {
            return members.to_vec();
        }
        panic!("members not defined")
    }

    fn get_node_id(&self, check_node: &Point) -> usize {
        if let Some(nodes) = &self.nodes {
            for node in nodes {
                if node == check_node {
                    return node.id;
                }
            }
        }

        panic!("Node does not exist");
    }

    fn get_node_id_special(&mut self, check_node: &Point) -> usize {
        if let Some(nodes) = &mut self.nodes {
            for node in nodes {
                let id = node.id;
                node.id = 0;
                if node == check_node {
                    node.id = id;
                    return node.id;
                }
                node.id = id;
            }
        }

        panic!("Node does not exist");
    }

    fn check_if_point_is_added(&mut self, point: &Point) -> bool {
        if let Some(nodes) = &mut self.nodes {
            // return nodes.contains(&point);
            for node in nodes {
                let id = node.id;
                node.id = 0;
                if node == point {
                    node.id = id;
                    return true;
                }
                node.id = id;
            }
        }
        false
    }

    fn add_node(&mut self, mut node: Point) -> usize {
        if self.check_if_point_is_added(&node) {
            return self.get_node_id_special(&node);
        }
        node.id = self.nodes_count;
        if let Some(nodes) = &mut self.nodes {
            nodes.push(node);
        } else {
            self.nodes = Some(vec![node]);
        }
        self.nodes_count += 1;
        self.nodes_count - 1
    }

    fn add_nodes(&mut self, nodes: Vec<Point>) -> Vec<usize> {
        let mut ids: Vec<usize> = vec![];
        for node in nodes {
            ids.push(self.add_node(node));
        }
        ids
    }

    fn add_member(&mut self, node1: Point, node2: Point) {
        if node1 == node2 {
            panic!("Two member Nodes cannot be same");
        }

        if self.check_if_point_is_added(&node1) && self.check_if_point_is_added(&node2) {
            let id1 = self.get_node_id_special(&node1);
            let id2 = self.get_node_id_special(&node2);

            if let Some(members) = &self.members {
                for member in members {
                    if member.contains(id1) && member.contains(id2) {
                        panic!(
                            "{} {}\n These two nodes are already added as members",
                            node1, node2
                        );
                    }
                }
            }
        }

        let ids = self.add_nodes(vec![node1, node2]);
        if let Some(members) = &mut self.members {
            members.push(Member::new(ids[0], ids[1], node1, node2));
            self.members_count+=1;
        } else {
            self.members = Some(vec![Member::new(ids[0], ids[1], node1, node2)]);
            self.members_count+=1;
        }
    }

    fn add_member_withparams(&mut self, node1: Point, node2: Point, e: f64, i: f64, a:f64) {
        if node1 == node2 {
            panic!("Two member Nodes cannot be same");
        }

        if self.check_if_point_is_added(&node1) && self.check_if_point_is_added(&node2) {
            let id1 = self.get_node_id_special(&node1);
            let id2 = self.get_node_id_special(&node2);

            if let Some(members) = &self.members {
                for member in members {
                    if member.contains(id1) && member.contains(id2) {
                        panic!(
                            "{} {}\n These two nodes are already added as members",
                            node1, node2
                        );
                    }
                }
            }
        }

        let ids = self.add_nodes(vec![node1, node2]);
        if let Some(members) = &mut self.members {
            members.push(Member::new_withparams(ids[0], ids[1], node1, node2, e, i, a));
        } else {
            self.members = Some(vec![Member::new_withparams(ids[0], ids[1], node1, node2, e, i, a)]);
        }
    }

    fn stiffness_matrix(&self) -> Matrix {
        let mut sm = Matrix::zeros(self.nodes_count*3);

        for (i, member) in  self.get_members().iter().enumerate(){
            let k = member.stiffness_matrix();
            let n = self.nodes_count;
            let x = k.columncount;
            let start = i*x/2 + n*i*x*x/4 as usize;
            let inc = n*x/2-x as usize;
            let mut factor: usize = 0;

            for j in 0..k.items.len() {
                sm.items[start+j+factor] += k.items[j];

                // println!("{}, {}", start+j+factor, j);

                if (j+1)%x == 0{
                    factor += inc
                }
            }

        }

        sm
    }

    fn add_roller_support(&mut self, point: Point, y_direction: bool) {
        let index = 3*self.get_node_id_special(&point);
        if let Some(vector) = &mut self.fixed_dofs {
            if y_direction {
                vector.push(index+1);
            }
            else {
                vector.push(index);
            }
        }

        else {
            if y_direction {
                self.fixed_dofs = Some(vec![index+1]);
            }
            else {
                self.fixed_dofs = Some(vec![index]);
            }
        }

    }

    fn add_hinge_support(&mut self, point: Point) {
        let index = 3*self.get_node_id_special(&point);
        if let Some(vector) = &mut self.fixed_dofs {
            vector.push(index);
            vector.push(index+1);
        }

        else {
            self.fixed_dofs = Some(vec![index, index+1]);
        }

    }

    fn add_fixed_support(&mut self, point: Point) {
        let index = 3*self.get_node_id_special(&point);
        if let Some(vector) = &mut self.fixed_dofs {
            vector.push(index);
            vector.push(index+1);
            vector.push(index+2);
        }

        else {
            self.fixed_dofs = Some(vec![index, index+1, index+2]);
        }

    }

    fn get_free_dofs(&self) -> Vec<usize> {
        let mut free_dofs: Vec<usize> = vec![];
        let fixed_dofs = self.get_fixed_dofs();

        for i in 0..3*self.nodes_count {
            if !fixed_dofs.contains(&i) {
                free_dofs.push(i)
            }
        }

        free_dofs
    }

    fn get_fixed_dofs(&self) -> Vec<usize> {
        if let Some(vector) = &self.fixed_dofs {
            return vector.to_vec();
        }

        panic!("Support not defined")
    }

    fn add_nodal_load(&mut self, x: f64, y: f64, m: f64, point: Point) {
        let index = 3*self.get_node_id_special(&point);
        if let Some(load_vector) = &mut self.load_vector {
            load_vector[index] += x;
            load_vector[index+1] += y;
            load_vector[index+2] += m;
        }
        else {
            let mut load_vector = vec![0.0; self.nodes_count*3];

            load_vector[index] = x;
            load_vector[index+1] = y;
            load_vector[index+2] = m;

            self.load_vector = Some(load_vector);
        }
    }

    fn add_varying_load(&mut self, begin: f64, end: f64, start_pt: Point, end_pt: Point) {
        let length = start_pt.distance(end_pt);
        let start_index = 3* self.get_node_id_special(&start_pt);
        let end_index = 3* self.get_node_id_special(&end_pt);

        let start_rxn = length * (7.0*begin + 3.0*end)/ 20.0;
        let end_rxn = length * (3.0*begin + 7.0*end)/ 20.0;

        let start_mom = -length*length * (-begin/20.0 -end/30.0);
        let end_mom = length*length * (-begin/30.0 -end/20.0);


        if let Some(load_vector) = &mut self.load_vector {
            
            load_vector[start_index+1] += start_rxn;
            load_vector[end_index+1] += end_rxn;

            load_vector[start_index+2] += start_mom;
            load_vector[end_index+2] += end_mom;
        }

        else {
            let mut load_vector = vec![0.0; self.nodes_count*3];

            load_vector[start_index+1] += start_rxn;
            load_vector[end_index+1] += end_rxn;

            load_vector[start_index+2] += start_mom;
            load_vector[end_index+2] += end_mom;

            self.load_vector = Some(load_vector);
        }
    }

    fn add_udl(&mut self, w: f64, start_pt: Point, end_pt: Point) {
        self.add_varying_load(w, w, start_pt, end_pt);
    }

    fn add_point_load_on_member(&mut self, load: f64, start_pt: Point, end_pt: Point, offset_from_start_pt: f64) {
        let length = start_pt.distance(end_pt);
        let start_index = 3* self.get_node_id_special(&start_pt);
        let end_index = 3* self.get_node_id_special(&end_pt);

        let b = offset_from_start_pt;
        let a = length-b;
        let start_rxn = load * a*a * (a+3.0*b)/(length*length*length);
        let end_rxn = load * b*b * (3.0*a+b)/(length*length*length);

        let start_mom = load*a*a*b/(length*length);
        let end_mom =  -load*a*b*b/(length*length);

        if let Some(load_vector) = &mut self.load_vector {
            
            load_vector[start_index+1] += start_rxn;
            load_vector[end_index+1] += end_rxn;

            load_vector[start_index+2] += start_mom;
            load_vector[end_index+2] += end_mom;
        }

        else {
            let mut load_vector = vec![0.0; self.nodes_count*3];

            load_vector[start_index+1] += start_rxn;
            load_vector[end_index+1] += end_rxn;

            load_vector[start_index+2] += start_mom;
            load_vector[end_index+2] += end_mom;

            self.load_vector = Some(load_vector);
        }
    }

    fn solve(&mut self) {
        if let Some(_) = &mut self.solution_load {
            self.solution_load = None;
            self.solve();
        }

        else {
            let load_vec = self.load_vector.clone();
            let mut vec_sol: Vec<f64> = vec![];
            if let Some(sol) = load_vec {
                vec_sol = sol.clone();
            }
            if let Some(_) = &mut self.load_vector {
                let gs = self.stiffness_matrix();
                let fv = vec_sol.clone();
                let free_dofs = self.get_free_dofs();
                let mut fixed_dofs = self.get_fixed_dofs();

                let mut rs = gs.clone();
                let mut rf = fv.clone();
                fixed_dofs.reverse();

                for i in fixed_dofs {
                    rs = rs.clone();
                    rs.remove_rowcol(i);
                    rf.remove(i);
                }


                let rd = rs.clone().inverse().vec_multipy(rf.clone());
                

                let mut count = 0;
                let mut dv = vec![0.0; self.nodes_count*3];
                for index in free_dofs {
                    dv[index] = rd[count];
                    count += 1;
                }

                let mut rv = gs.clone().vec_multipy(dv.clone());

                for i in 0..rv.len() {
                    rv[i] -= fv[i];
                }

                self.solution_displacement = Some(dv.clone());
                self.solution_load = Some(rv.clone());

                println!("--------------------------------");
                println!("Solved Forces   |  Displacements");
                println!("--------------------------------");
                for (i, (f, d)) in (rv.iter().zip(dv)).into_iter().enumerate() {
                    if (i) % 3 == 0 {
                        println!("H{} = {}\t|\t\u{0394}x = {}", i/3+1, format!("{:.3}", f), format!("{:.4}", d));
                    }
                    else if (i-1) % 3 == 0 {
                        println!("V{} = {}\t|\t\u{0394}y = {}", i/3+1, format!("{:.3}", f), format!("{:.4}", d));
                    }
                    else { 
                        println!("M{} = {}\t|\t\u{03B8}  = {}", i/3+1, format!("{:.3}", f), format!("{:.4}", d));
                    }

                    if (i+1)%3 == 0 {
                        println!("--------------------------------");
                    }
                }
            }
        }

    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_count() {
        let a = Point::new(1.0, 2.0);
        let b = Point::new(-1.0, 4.0);
        let c = Point::new(9.3, 3.3);
        let d = Point::new(0.3, 0.2);
        let e = Point::new(9.0, 2.1);
        let fa = Point::new(8.3, 3.4);
        let mut f = Frame::new();

        f.add_member(a, b);
        f.add_member(b, c);
        f.add_member(c, d);
        f.add_member(d, e);
        f.add_member(fa, a);
        f.add_member(e, fa);
        assert_eq!(f.nodes_count, 6)
    }

    #[test]
    fn distance_test() {
        let p1 = Point::new_from_int(0, 3);
        let p2 = Point::new_from_int(4, 0);
        assert_eq!(p1.distance(p2), 5.0)
    }

    #[test]
    fn stiff_test() {
        let a = Point::new_from_int(0, 0);
        let b = Point::new_from_int(8, 0);
        let c = Point::new_from_int(14, 0);

        let mut f = Frame::new();

        f.add_member_withparams(a, b, 20e5, 5e-3, 1.0);
        f.add_member_withparams(b, c, 20e5, 5e-3, 1.0);
        let k = f.stiffness_matrix();

        let  stiff = vec![250000.0, 0.0, 0.0, -250000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 234.375, 937.5, 0.0, -234.375, 937.5, 0.0, 0.0, 0.0, 0.0, 937.5, 5000.0, 0.0, -937.5, 2500.0, 0.0, 0.0, 0.0, -250000.0, 0.0, 0.0, 583333.3333333333, 0.0, 0.0, -333333.3333333333, 0.0, 0.0, 0.0, -234.375, -937.5, 0.0, 789.9305555555555, 729.1666666666667, 0.0, -555.5555555555555, 1666.6666666666667, 0.0, 937.5, 2500.0, 0.0, 729.1666666666667, 11666.666666666668, 0.0, -1666.6666666666667, 3333.3333333333335, 0.0, 0.0, 0.0, -333333.3333333333, 0.0, 0.0, 333333.3333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -555.5555555555555, -1666.6666666666667, 0.0, 555.5555555555555, -1666.6666666666667, 0.0, 0.0, 0.0, 0.0, 1666.6666666666667, 3333.3333333333335, 0.0, -1666.6666666666667, 6666.666666666667];

        assert_eq!(k.items, stiff);
        
    }

    #[test]
    fn dof_test() {

        let a = Point::new_from_int(0, 0);
        let b = Point::new_from_int(8, 0);
        let c = Point::new_from_int(14, 0);

        let mut f = Frame::new();

        f.add_member_withparams(a, b, 20e5, 5e-3, 1.0);
        f.add_member_withparams(b, c, 20e5, 5e-3, 1.0);

        f.add_fixed_support(a);
        f.add_fixed_support(c);

        assert_eq!(f.get_fixed_dofs(), vec![0, 1, 2, 6, 7, 8]);
        assert_eq!(f.get_free_dofs(), vec![3, 4, 5]);
    }

    #[test]
    fn load_test() {
        let a = Point::new_from_int(0, 0);
        let b = Point::new_from_int(8, 0);
        let c = Point::new_from_int(14, 0);

        let mut f = Frame::new();

        f.add_member_withparams(a, b, 20e5, 5e-3, 1.0);
        f.add_member_withparams(b, c, 20e5, 5e-3, 1.0);

        f.add_fixed_support(a);
        f.add_roller_support(b, true);
        f.add_roller_support(c, true);

        f.add_point_load_on_member(-30.0, a, b, 3.0);
        f.add_nodal_load(0.0, 0.0, -40.0, b);
        f.add_udl(-20.0, b, c);

        assert_eq!(f.load_vector, Some(vec![0.0, -20.5078125, -35.15625, 0.0, -69.4921875, -78.90625, 0.0, -60.0, 59.99999999999999]));
    }
}

fn main() {

    let file_path = "./src/data";
    let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    let mut points: Vec<PData> = vec![];
    let mut members: Vec<MData> = vec![];
    let mut nodalloads: Vec<NodalLoad> = vec![];
    let mut pointloads: Vec<PointLoad> = vec![];
    let mut udls: Vec<UDL> = vec![];
    let mut uvls: Vec<UVL> = vec![];
    let mut supports: Vec<SData> = vec![];

    let mut point_read_start_flag = false;
    let mut members_read_start_flag = false;
    let mut nodal_load_read_start_flag = false;
    let mut point_load_read_start_flag = false;
    let mut udl_read_start_flag = false;
    let mut uvl_read_start_flag = false;
    let mut support_read_start_flag = false;

    // parsing the data from file_path
    for line in contents.split("\n") {
        if line.contains("//") {continue;}

        if line.contains("[Supports]") {
            point_read_start_flag = false;
            members_read_start_flag = false;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = false;
            udl_read_start_flag = false;
            uvl_read_start_flag = false;
            support_read_start_flag = true;
            continue;
        }

        if line.contains("[UVL]") {
            point_read_start_flag = false;
            members_read_start_flag = false;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = false;
            udl_read_start_flag = false;
            uvl_read_start_flag = true;
            support_read_start_flag = false;
            continue;
        }

        if line.contains("[UDL]") {
            point_read_start_flag = false;
            members_read_start_flag = false;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = false;
            udl_read_start_flag = true;
            uvl_read_start_flag = false;
            support_read_start_flag = false;
            continue;
        }

        if line.contains("[PointLoads]") {
            point_read_start_flag = false;
            members_read_start_flag = false;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = true;
            udl_read_start_flag = false;
            uvl_read_start_flag = false;
            support_read_start_flag = false;
            continue;
        }

        if line.contains("[NodalLoads]") {
            point_read_start_flag = false;
            members_read_start_flag = false;
            nodal_load_read_start_flag = true;
            point_load_read_start_flag = false;
            udl_read_start_flag = false;
            uvl_read_start_flag = false;
            support_read_start_flag = false;
            continue;
        }

        if line.contains("[Members]") {
            point_read_start_flag = false;
            members_read_start_flag = true;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = false;
            udl_read_start_flag = false;
            uvl_read_start_flag = false;
            support_read_start_flag = false;
            continue;
        }

        if line.contains("[Points]") {
            point_read_start_flag = true;
            members_read_start_flag = false;
            nodal_load_read_start_flag = false;
            point_load_read_start_flag = false;
            udl_read_start_flag = false;
            uvl_read_start_flag = false;
            support_read_start_flag = false;
            continue;
        }

        if support_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = SData::new("".to_string(), 0, true);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.typ = pt.trim().to_string();
                }
                if i == 1 {
                    point.index = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.y_orient = pt.trim().parse().unwrap();
                }

            }
            supports.push(point);
        }

        if uvl_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = UVL::new(0.0, 0.0, 0, 0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.start_val = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.end_val = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.start = pt.trim().parse().unwrap();
                }
                if i == 3 {
                    point.end = pt.trim().parse().unwrap();
                }
            }
            uvls.push(point);
        }

        if udl_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = UDL::new(0.0, 0, 0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.val = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.start = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.end = pt.trim().parse().unwrap();
                }
            }
            udls.push(point);
        }

        if point_load_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = PointLoad::new(0.0, 0, 0, 0.0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.val = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.start = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.end = pt.trim().parse().unwrap();
                }
                if i == 3 {
                    point.offset = pt.trim().parse().unwrap();
                }
            }
            pointloads.push(point);
        }

        if nodal_load_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = NodalLoad::new(0.0, 0.0, 0.0, 0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.x = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.y = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.m = pt.trim().parse().unwrap();
                }
                if i == 3 {
                    point.index = pt.trim().parse().unwrap();
                }
            }
            nodalloads.push(point);
        }

        if members_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = MData::new(0, 0, 1.0, 1.0, 1.0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.x = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.y = pt.trim().parse().unwrap();
                }
                if i == 2 {
                    point.e = pt.trim().parse().unwrap();
                }
                if i == 3 {
                    point.i = pt.trim().parse().unwrap();
                }
                if i == 4 {
                    point.a = pt.trim().parse().unwrap();
                }
            }
            members.push(point);
        }

        if point_read_start_flag {
            let line = line.trim();
            if line == "" {
                continue;
            }

            let mut point = PData::new(0.0, 0.0);
            for (i, pt) in line.split(",").into_iter().enumerate() {
                if i == 0 {
                    point.x = pt.trim().parse().unwrap();
                }
                if i == 1 {
                    point.y = pt.trim().parse().unwrap();
                }
            }
            points.push(point);
        }

    }

    let mut points_conv: Vec<Point> = vec![];

    for point in points {
        points_conv.push(Point::new(point.x, point.y));
    }

    let mut frame = Frame::new();

    for m in members {
        frame.add_member_withparams(
            points_conv[m.x], points_conv[m.y], 
            m.e, m.i, m.a
            );
    }

    for s in supports {
        if s.typ == "roller" {
            frame.add_roller_support(points_conv[s.index], s.y_orient);
        }
        else if s.typ == "hinge" {
            frame.add_hinge_support(points_conv[s.index]);
        }
        else if s.typ == "fixed" {
            frame.add_fixed_support(points_conv[s.index]);
        }
    }

    for nl in nodalloads {
        frame.add_nodal_load(nl.x, nl.y, nl.m, points_conv[nl.index])
    }
    for pl in pointloads{
        frame.add_point_load_on_member(pl.val,
            points_conv[pl.start],
            points_conv[pl.end],
            pl.offset)
    }
    for u in udls {
        frame.add_udl(
            u.val,
            points_conv[u.start],
            points_conv[u.end],
        )
    }
    for v in uvls {
        frame.add_varying_load(
            v.start_val,
            v.end_val,
            points_conv[v.start],
            points_conv[v.end],
        )
    }

    frame.solve();
}
