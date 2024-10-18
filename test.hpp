#ifndef TEST_HPP
#define TEST_HPP

using std::cout;

void test1() {

    // testing line intersect
    cout << "testing line intersect" <<endl;
    Vec<double> line0{1.0, -1.0, 3.0};
    Vec<double> line1{-2.0, 1.0, 4.0};
    Vec<double> line30{0.0, 1.0, -1.0};
    Vec<double> line31{1.0, 0.0, -0.5};
    //cout << line_line_intersection(line0, line1) << endl;
    cout << line_line_intersection(line30, line31) << endl;
    cout << "done testing line intersect" << endl;
    // testing segment to line
    cout << "testing segment to line" << endl;
    Vec<double> p0{2.0, 5.0};
    Vec<double> p1{10.0, 25.0};
    cout << get_line_from_segment(p0, p1) << endl;
    cout << "done testing segment to line" << endl;
    // testing intersection between segment and line
    cout << "testing segment to line intersect" << endl;
    cout << "testing segment to line intersect" << endl;
    Vec<double> line2{1.0, 0.0, -5.0};
    Vec<double> line20{0.0, 1.0, -1.0};
    Vec<double> p100{0.5, -0.5};
    Vec<double> p101{0.5, 2.0};
    cout << "  double check segment to line " << get_line_from_segment(p100, p101) << endl;
    cout << " line20 " << line20 << endl;
    cout << " p100 " << p100 << endl << " p101 " << p101 << endl;
    //cout << line_line_segment_intersection(line2, p0, p1) << endl;
    cout << " " << line_line_segment_intersection(line20, p100, p101) << endl;
    cout << "done testing segment to line intersect" << endl;
    //testing point in region
    cout << "point in region" << endl;
    Vec<double> line3{1.0, 0.0, -5.0};
    Vec<double> p2{-10.0, 1.0};
    Vec<double> p3{100.0, -3.0};
    Vec<double> p4{1.0, 100.0};
    Vec<double> p5{1.0, 100.0};
    // testing slicing again
    cout << "testing slicing" << endl;
    Vec<double> p6{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; 
    Vec<double> p7 = p6.slice(0, 2);
    Vec<double> p8 = p6.slice(1, 3);
    Vec<double> p9 = p6.slice(2, 4);
    //Vec<double> p10 = p6.slice(4, 3);
    cout << "p6 " << p6 << endl << p7 << endl << p8 << endl << p9 << endl;
    cout << "done testing slicing" << endl;
    // testing clipping
    using std::vector;
    cout << "testing clip polygon function" << endl;
    vector<Vec<double>> input_vertices = { Vec{0.5, 2.0}, Vec{0.5, -0.5}, Vec{-0.5,-0.5}, Vec{-0.5, 0.5} };
    vector<Vec<double>> edges = { Vec{1.0, 0.0, -1.0}, Vec{0.0, 1.0, 1.0}, Vec{1.0, 0.0, 1.0}, Vec{0.0, 1.0, -1.0} };
    Vec<double> pt_inside = {0.0, 0.0};
    vector<Vec<double>> output_vertices(8);
    clip_polygon(input_vertices, edges, pt_inside, output_vertices);
    cout << " printing input polygon vertices" << endl;
    for (auto& input_vertices : input_vertices) {
        cout << " " << input_vertices << endl;
    }
    cout << " printing output polygon vertices" << endl;
    for (auto& output_vertex : output_vertices) {
        cout << " " << output_vertex << endl;
    }
    cout << " done printing clip polygon vertices" << endl;
    cout << point_in_region(line1, p2, p3) << endl;
    cout << point_in_region(line1, p2, p4) << endl;
    cout << point_in_region(line3, p2, p3) << endl;
    cout << point_in_region(line3, p2, p4) << endl;
    cout << "done point in region" << endl;
    // testing object compare
    cout << "p4 == p5 " << (p4 == p5) << endl;
    cout << "p4 == p3 " << (p4 == p3) << endl;
    // testing line intersect plane
    cout << "Testing line intersect plane" << endl;
    Vec line_a{2.0, -1.0, 2.0, 1.0, 4.0, 1.0};
    Vec plane_a{1.0, -2.0, 1.0, -12.0};
    cout << line_plane_intersect(line_a, plane_a) << endl;
    cout << "done testing line intersect plane" << endl;
    // testing point in region
    cout << "testing point in region3d" << endl;
    Vec pa{-1.0, 5.0, 10.0};
    Vec pb{50.0,3.0,-5.0};
    Vec pc{-100.0, -10.0, -20.0};
    Vec pd{75.0, -100.0, 1050.0};
    Vec plane_c{1.0, 0.0, 0.0, -3.0};
    cout << point_in_region3d(plane_c, pa, pb) << endl; // false
    cout << point_in_region3d(plane_c, pa, pc) << endl; // true
    cout << point_in_region3d(plane_c, pa, pd) << endl; // false
    cout << "done testing point in region3d" << endl;
    //T dist_from_plane(Vec<T>& plane, Vec<T>& point)
    // testing distance from plane
    cout << "testing distance from plane" << endl;
    Vec pf{8.0, 0.0, 0.0};  
    Vec pg{-2.0, 5.0, 0.0};  
    cout << dist_from_plane(plane_c, pf) << endl;
    cout << dist_from_plane(plane_c, pg) << endl;
    cout << "done testing distance from plane" << endl;
    // testing point on line 3d
    cout << "testing pt on line3d\n";
    Vec linea3d{1.0,2.0,1.0,2.0,4.0,0.0};
    Vec pa3d{3.0, 6.0,0.0};
    cout << "on line? " << point_on_line3d(pa3d, linea3d) << endl;
    //testing point on line segment 3d
    cout << "testing pt on line segment 3d\n";
    Vec pseg0{2.0,2.0,2.0};
    Vec pseg1{6.0,2.0,2.0};
    Vec pseg2{4.0,2.0,2.0};
    Vec pseg3{4.0,3.0,2.0};
    cout << "  should be 1: \n" << point_on_line_segment3d(pseg2, pseg0, pseg1) << endl;
    cout << "  should be 0: \n" << point_on_line_segment3d(pseg3, pseg0, pseg1) << endl;
   
        // new segment
    Vec pseg4{ 27.0, 36.0,-6.0};
    Vec pseg5{-17.0, -52.0, 27.0}; // (t = -6)
        // test points on new segment
    Vec pseg6{-13.0,-44.0,24.0};  // t = -5
    Vec pseg7{-20.0, -55.0, 24.0}; // should not be on line
    cout << "  should be 1: \n" << point_on_line_segment3d(pseg6, pseg4, pseg5) << endl;
    cout << "  should be 0: \n" << point_on_line_segment3d(pseg7, pseg4, pseg5) << endl;
    // test plane line segment intersect
    cout << "testing plane line segment intersect\n";
    cout << "testing plane line segment intersect\n";
    cout << "testing plane line segment intersect\n";
    Vec plane0{1.0,0.0,0.0,-4.0};
    Vec pseg10{-10.0,0.0,0.0};
    Vec pseg11{20.0,0.0,0.0};
    Vec pseg12{-100.0,0.0,0.0};
    cout << "  4th element should be nonzero" << plane_line_segment_intersect(plane0, pseg10, pseg11) << endl;
    cout << "  4th element should be zero" << plane_line_segment_intersect(plane0, pseg10, pseg12) << endl;
        // new segment 
    Vec plane1{1.0, 3.0, -2.0, 4.0};
        // this segment passes thru plane
    Vec pseg13{126.0/41.0, -109.0/41.0, -50.0/41.0 }; 
    Vec pseg14{0.0/41.0, -109/41.0, -50.0/41.0 };
        // this segment does not pass thru the plane
    Vec pseg15{11.5, 14.5, 7.0};
    Vec pseg16{8.5, 5.5, 13.0};
    cout << "  4th element should be nonzero" << plane_line_segment_intersect(plane1, pseg13, pseg14) << endl;
    cout << "  4th element should be zero" << plane_line_segment_intersect(plane1, pseg15, pseg16) << endl;
        // new segment and plane
    Vec pseg17{ 0.0402188, -0.0402188, 1.00001 };
    Vec pseg18{ -1.15741,  -1.15741,   1.00001 };
    Vec plane2{1.0, 0.0, 0.0, 1.0};
    cout << "  hi there " << plane_line_segment_intersect(plane2, pseg17, pseg18) << endl;


    cout << "tests finished" << endl;
}

// older tests, mostly matrix and vector operations
void test0() {
    //float temp = bisection( f, 1.0, 2.0);
    //cout << "bisection = " << temp << endl;
    //temp = false_position( f, 1.0, 2.0);
    //cout << "false_position = " << temp << endl;

    double c = 10.0;
    Vec<double> a = {1.0, 2.0, 3.0, 4.0 };
    //Vec<double> a(4); 
    //a.data[0] = 1.0;
    //a.data[1] = 2.0;
    //a.data[2] = 3.0;
    //a.data[3] = 4.0;
    cout << " slicing a: index 1, count 2 : " << a.slice(1,2) << endl;
    Vec<double> b(4);
    b.data[0] = 10.0;
    b.data[1] = 20.0;
    b.data[2] = 30.0;
    b.data[3] = 40.0;
    cout << a << endl;
    cout << b << endl;
    cout << c * b  << endl;
    cout << a.mag() << endl;
    cout << a.unit() << endl;
    cout << a.dot(b) << endl;
    Vec<double> a1(3); 
    a1.data[0] = -20.0;
    a1.data[1] = 2.0;
    a1.data[2] = 3.0; // change to 4.0
    Vec<double> b1(3);
    b1.data[0] = 10.0;
    b1.data[1] = 20.0;
    b1.data[2] = 30.0;
    cout << "a1 cross b1 " << a1.cross(b1) << endl; 
    // test matrix 
    Mat<double> m(4, 4);
    m.set_row(0, a); //m.rows[0] = a;
    m.set_row(1, b); //m.rows[1] = b;
    m.set_row(2, a + b); //m.rows[2] = a + b;
    cout << m << endl; 
    //print_matrix(m);

    m.set(1, 3, 780.2);
    m.set(3, 2, 80.0);
    cout << m << endl; 
    cout << "scaling matrix by 3.0" << endl;
    cout << 3.0 * m << endl;
    Mat<double> m1(4, 4);
    m1.set_row(0, a); //m1.rows[0] = a;
    m1.set_row(1, b); //m1.rows[1] = b;
    m1.set_row(2, a + b); //m1.rows[2] = a + b;
    m1.set(1, 2, -28.0);
    m1.set(0, 1, -40.0);
    m1.set(3, 1, 15.0);
    cout << "m" << endl << m << endl;
    cout << "m1" << endl << m1 << endl;
    cout << "testing matrix multiply" <<endl;
    Mat<double> mat_new(m * m1); // copy const
    cout << mat_new << endl;
    cout << "testing matrix transpose of mat_new" << endl;
    cout << mat_new.transpose() << endl; 
    cout << "sub matrix m" << endl;
    cout << m.sub_matrix(1, 2) << endl;
    m.set(0,3, -48.0);
    m.set(3,0, 1.0);
    m.set(1,1, -7.0);
    m.set(3,1, -18.0);
    m.set(3,3, 17.0);
    cout << "reprinting edited matrix m" << endl;
    cout << m << endl;
    cout << "printing determinant of m" << endl;
    cout << m.det() << endl;
    cout << "inverse of m" << endl;
    cout << m.inv() << endl;
    Mat<double> m2(4, 4);
    for (int i = 0; i < 16; ++i) {
        m2.data[i] = 2*i;
    }
    m2.data[0] = 3.0;
    Vec<double> v2(4);
    for(int i = 0; i < 4; ++i) {
        v2.data[i] = 2*i + 3;
    }
    cout << "m2\n" << m2 << "v2\n" << v2 << endl;
    cout << "vector * matrix" << endl;
    cout << v2 * m2 << endl;
    
}
#endif
