// cd project3
// cd part4
// g++ -std=c++11 -o l034  -Wall l034.cpp
// ./l034

#include <iostream>
#include <fstream>
#include <iomanip>  // set precision
#include <climits>  // int_max
#include <list>     // list
#include <math.h>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>    // std::sort
#include <map>
#include <unordered_map>
#include <functional>

const int size = 800;
std::ofstream results;
double duration1;
double duration2;
double duration3;
double duration4;

// point class
// private, getters, setters, constructors
class Point {
public:
    void setX(double x);
    double getX(void);
    void setY(double y);
    double getY(void);
    bool diff(Point p2);
    bool same(Point p2);
    //int Compare(Point d);
    static bool sortFunct(Point i, Point j) { return (j.getX() > i.getX()); }
    static std::vector<Point> getSorted(std::vector<Point> p) {
        std::sort(p.begin(), p.end(), sortFunct);
        return p;
    }
    static bool sortFunctY(Point i, Point j) { return (j.getY() > i.getY()); }
    static std::vector<Point> getSortedY(std::vector<Point> p) {
        std::sort(p.begin(), p.end(), sortFunctY);
        return p;
    }
    void print();
    Point(double x, double y);  // constructor
    Point(void);
private:
    double x;
    double y;
};
Point::Point(void) {
    x = 0;
    y = 0;
}
Point::Point(double xval, double yval) {
    x = xval;
    y = yval;
}
void Point::setX(double xval) {
    x = xval;
}
double Point::getX(void) {
    return x;
}
void Point::setY(double yval) {
    y = yval;
}
double Point::getY(void) {
    return y;
}
bool Point::diff(Point p2) {
    if (getX() != p2.getX() && getY() != p2.getY())
        return true;
    return false;
}
bool Point::same(Point p2) {
    if (getX() == p2.getX() && getY() == p2.getY())
        return true;
    return false;
}
void Point::print() {
    std::cout << std::fixed;
    std::cout << std::setprecision(23) << "(" << getX() << " , " << getY() << ")";
}

class Distance {
public:
    void setdist(double d);
    double getdist(void);
    void setp1(Point p);
    Point getp1(void);
    void setp2(Point p);
    Point getp2(void);
    Distance compare(Distance d1, Distance d2);
    void print();
    Distance(double dist, Point p1, Point p2);  // constructors
    Distance(void);
private:
    double dist;
    Point p1;
    Point p2;
};
Distance::Distance(void) {
    dist = 0;
    p1 = Point();
    p2 = Point();
}
Distance::Distance(double d, Point p1val, Point p2val) {
    dist = d;
    p1 = p1val;
    p2 = p2val;
}
void Distance::setdist(double d) {
    dist = d;
}
double Distance::getdist(void) {
    return dist;
}
void Distance::setp1(Point p) {
    p1 = p;
}
Point Distance::getp1(void) {
    return p1;
}
void Distance::setp2(Point p) {
    p2 = p;
}
Point Distance::getp2(void) {
    return p2;
}
Distance Distance::compare(Distance d1, Distance d2) {
    if (d1.getdist() < d2.getdist())
        return d1;
    return d2;
}
void Distance::print() {
    std::cout << std::fixed;
    std::cout << std::setprecision(23) << dist << "\t";
    getp1().print();
    std::cout << "\t";
    getp2().print();
}

class Box {
public:
    void setX(unsigned long long int x);
    unsigned long long int getX(void)const;
    void setY(unsigned long long int y);
    unsigned long long int getY(void)const;
    bool operator==(const Box& b) const;
    void print();
    Box(unsigned long long int x, unsigned long long int y);  // constructor
    Box(void);
private:
    unsigned long long int x;
    unsigned long long int y;
};
Box::Box(void) {
    x = 0;
    y = 0;
}
Box::Box(unsigned long long int xval, unsigned long long int yval) {
    x = xval;
    y = yval;
}
void Box::setX(unsigned long long int xval) {
    x = xval;
}
unsigned long long int Box::getX(void) const {
    return x;
}
void Box::setY(unsigned long long int yval) {
    y = yval;
}
unsigned long long int Box::getY(void) const {
    return y;
}
bool Box::operator==(const Box& b) const {
    return x == b.x && y == b.y;
}
void Box::print() {
    std::cout << std::fixed;
    std::cout << std::setprecision(23) << "(" << getX() << " , " << getY() << ")";
}

class MyHashFunction {
public:
    // We use predefined hash functions of strings
    // and define our hash function as XOR of the
    // hash values.
    size_t operator()(const Box& b) const {
        return (std::hash<unsigned long long int>()(b.getX())) ^ (std::hash<unsigned long long int>()(b.getY()));
    }
};

void illuminate(int** arr, int y, int x, int color) {
    if (color == 0)
        arr[x][y] = 0;
    else if (color == 1)
        arr[x][y] = 1;
    else if (color == 2)
        arr[x][y] = 2;
    else if (color == 3)
        arr[x][y] = 3;
}

void notoutofbounds(int** array, int x, int y, int color) {
    if (x >= 0 && x < size && y >= 0 && y < size)
        //array[x][y] = 0;
        illuminate(array, x, y, color);
}

double distance(double x1, double y1, double x2, double y2) {
    double t1 = pow(x1 - x2, 2);
    double t2 = pow(y1 - y2, 2);
    return (double)(sqrt(t1 + t2));
}

void circle(int** array, int xcenter, int ycenter, int R, int color) {
    int y = R;
    int xmax = (int)(R * 0.70710678);
    int ysquare = y * y;
    int ty = (2 * y) - 1;
    int y2_new = ysquare;

    //std::cout << y << " " << xmax << " " << ysquare << " " << ty << " " << y2_new << "\n";

    for (int x = 0; x <= xmax + 1; x++) {
        if ((ysquare - y2_new) >= ty) {
            ysquare -= ty;
            y -= 1;
            ty -= 2;
        }
        notoutofbounds(array, y + xcenter, x + ycenter, color);
        notoutofbounds(array, -y + xcenter, x + ycenter, color);
        notoutofbounds(array, y + xcenter, -x + ycenter, color);
        notoutofbounds(array, -y + xcenter, -x + ycenter, color);

        notoutofbounds(array, x + xcenter, y + ycenter, color);
        notoutofbounds(array, -x + xcenter, y + ycenter, color);
        notoutofbounds(array, x + xcenter, -y + ycenter, color);
        notoutofbounds(array, -x + xcenter, -y + ycenter, color);

        y2_new -= (2 * x) - 3;
    }
}

Distance part1() {
    std::ofstream points;
    points.open("points.txt");
    std::list <Point*> pointlist;

    // make array
    int** array;
    array = new int* [size];

    // fill with 1s
    for (int row = 0; row < size; ++row)
        array[row] = new int[size];
    for (int row = 0; row < size; ++row)
        for (int col = 0; col < size; ++col)
            array[row][col] = 1;

    srand(time(NULL));
    for (int i = 0; i <= 999; i++) {
        double x = double(rand()) / RAND_MAX;
        double y = double(rand()) / RAND_MAX;
        pointlist.push_back(new Point(x, y));
        circle(array, (int)(x * size), (int)(y * size), 3, 0);
        points << std::fixed;
        points << std::setprecision(23) << x << "  " << y << "\n";
    }
    points.close();

    // make ppm
    std::ofstream ppm;
    ppm.open("points1.ppm");
    ppm << "P3 " << size << " " << size << " 255\n";

    std::clock_t start = std::clock();

    // iterator
    double d = INT_MAX;
    Point dp1;
    Point dp2;
    Point temp1;
    Point temp2;
    std::list <Point*> ::iterator it1;
    std::list <Point*> ::iterator it2;
    for (it1 = pointlist.begin(); it1 != pointlist.end(); it1++) {
        for (it2 = it1++; it2 != pointlist.end(); it2++) {
            temp1 = **it1;
            temp2 = **it2;
            if (temp1.diff(temp2)) {
                double dist = distance(temp1.getX(), temp1.getY(), temp2.getX(), temp2.getY());
                if (dist < d) {
                    d = dist;
                    dp1 = temp1;
                    dp2 = temp2;
                }
            }
        }
    }

    double duration = ((double)std::clock() - (double)start) / (double)CLOCKS_PER_SEC;
    duration1 = duration;

    circle(array, (int)(dp1.getX() * size), (int)(dp1.getY() * size), 2, 2);
    circle(array, (int)(dp1.getX() * size), (int)(dp1.getY() * size), 3, 2);
    circle(array, (int)(dp2.getX() * size), (int)(dp2.getY() * size), 2, 2);
    circle(array, (int)(dp2.getX() * size), (int)(dp2.getY() * size), 3, 2);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (array[i][j] == 0)
                ppm << array[i][j] << " " << array[i][j] << " " << array[i][j] << ' ';
            else if (array[i][j] == 2)
                ppm << 255 << " " << 0 << " " << 0 << " ";
            else
                ppm << 255 << " " << 255 << " " << 255 << ' ';
        }
        ppm << std::endl;
    }
    ppm << "\n";

    ppm.close();

    Distance dist = Distance(d, Point(dp1.getX(), dp1.getY()), Point(dp2.getX(), dp2.getY()));

    return dist;
}

Distance stripClosest(std::vector<Point> pointlist) {
    double d = INT_MAX;
    Point dp1;
    Point dp2;
    Point temp1;
    Point temp2;
    for (size_t i = 0; i < pointlist.size(); i++) {
        for (size_t j = 0; j < pointlist.size(); j++) {
            temp1 = pointlist[i];
            temp2 = pointlist[j];
            if (temp1.diff(temp2)) {
                double dist = distance(temp1.getX(), temp1.getY(), temp2.getX(), temp2.getY());
                if (dist < d) {
                    d = dist;
                    dp1 = temp1;
                    dp2 = temp2;
                }
            }
        }
    }
    return Distance(d, dp1, dp2);
}

Distance recursive(std::vector<Point> v) {
    std::vector<Point> left;
    std::vector<Point> right;

    for (int i = 0; i != (int)v.size() / 2; i++)
        left.push_back(v[i]);
    for (int i = v.size() / 2; i < (int)v.size(); i++)
        right.push_back(v[i]);

    if (left.size() == 2 || left.size() == 3) { // base case left strip
        double dl;
        double d1 = distance(left[0].getX(), left[0].getY(), left[1].getX(), left[1].getY());
        if (left.size() == 3) {
            double d2 = distance(left[1].getX(), left[1].getY(), left[2].getX(), left[2].getY());
            double d3 = distance(left[0].getX(), left[0].getY(), left[2].getX(), left[2].getY());
            dl = std::min(d3, std::min(d1, d2));
            if (dl == d2) {
                return Distance(dl, left[1], left[2]);
            }
            if (dl == d3) {
                return Distance(dl, left[0], left[2]);
            }
        }
        dl = d1;
        return Distance(dl, left[0], left[1]);
    }
    if (right.size() == 2 || right.size() == 3) { // base case right strip
        double dr;
        double d1 = distance(right[0].getX(), right[0].getY(), right[1].getX(), right[1].getY());
        if (right.size() == 3) {
            double d2 = distance(right[1].getX(), right[1].getY(), right[2].getX(), right[2].getY());
            double d3 = distance(right[0].getX(), right[0].getY(), right[2].getX(), right[2].getY());
            dr = std::min(d3, std::min(d1, d2));
            if (dr == d2) {
                return Distance(dr, right[1], right[2]);
            }
            if (dr == d3) {
                return Distance(dr, right[0], right[2]);
            }
        }
        dr = d1;
        return Distance(dr, right[0], right[1]);
    }

    // Find the middle point
    int mid = v.size() / 2;
    Point midPoint = v[mid];

    // Consider the vertical line passing through the middle point calculate
    // the smallest distance dl on left of middle point and dr on right side
    Distance dl = recursive(left);
    Distance dr = recursive(right);

    // Find the smaller of two distances
    Distance d = dl.compare(dl, dr);

    // Build an array strip[] that contains points close (closer than d)
    // to the line passing through the middle point
    std::vector<Point> strip;
    for (size_t i = 0; i < v.size(); i++)
        if (abs(v[i].getX() - midPoint.getX()) < d.getdist())
            strip.push_back(v[i]);

    // Find the closest points in strip. Return the minimum of d and closest distance is strip[]
    return d.compare(d, stripClosest(strip));
}

Distance part2() {
    //std::cout << "@1\n";
    std::vector<Point> v = {};

    //pulling the numbers from the file
    std::vector<std::string> vec = {};
    std::string line;
    std::ifstream points("points.txt");
    //std::cout << "@2\n";

    if (points.is_open()) {
        //std::cout << "@3\n";
        while (getline(points, line)) {
            std::stringstream ss(line);
            while (getline(ss, line, ' ')) {
                if (line.length() > 2)
                    vec.push_back(line);
            }
        }
        //std::cout << "@4\n";
        points.close();
    }
    else
        std::cout << "Unable to open file";
    //std::cout << "@5\n";

    //std::cout << "@5.5\n";
    std::vector<std::string>::iterator it;
    //std::cout << "@6\n";
    for (it = vec.begin(); it != vec.end(); it++) {
        std::string s1 = *it;
        //std::cout << *it << " first value\n";
        it++;
        //std::cout << *it << " second value\n";
        std::string s2 = *it;
        v.push_back(Point(std::stod(s1), std::stod(s2)));
        //v.back().print();
    }

    //std::cout << "@7\n";

    std::clock_t start = std::clock();

    /* Your algorithm m */
    v = Point::getSorted(v);

    Distance d = recursive(v);
    double duration = ((double)std::clock() - (double)start) / (double)CLOCKS_PER_SEC;
    duration2 = duration;

    return d;
}

Distance recursive3(std::vector<Point> v) {
    std::vector<Point> left;
    std::vector<Point> right;

    for (int i = 0; i != (int)v.size() / 2; i++)
        left.push_back(v[i]);
    for (int i = v.size() / 2; i < (int)v.size(); i++)
        right.push_back(v[i]);

    if (left.size() == 2 || left.size() == 3) { // base case left strip
        double dl;
        double d1 = distance(left[0].getX(), left[0].getY(), left[1].getX(), left[1].getY());
        if (left.size() == 3) {
            double d2 = distance(left[1].getX(), left[1].getY(), left[2].getX(), left[2].getY());
            double d3 = distance(left[0].getX(), left[0].getY(), left[2].getX(), left[2].getY());
            dl = std::min(d3, std::min(d1, d2));
            if (dl == d2) {
                return Distance(dl, left[1], left[2]);
            }
            if (dl == d3) {
                return Distance(dl, left[0], left[2]);
            }
        }
        dl = d1;
        return Distance(dl, left[0], left[1]);
    }
    if (right.size() == 2 || right.size() == 3) { // base case right strip
        double dr;
        double d1 = distance(right[0].getX(), right[0].getY(), right[1].getX(), right[1].getY());
        if (right.size() == 3) {
            double d2 = distance(right[1].getX(), right[1].getY(), right[2].getX(), right[2].getY());
            double d3 = distance(right[0].getX(), right[0].getY(), right[2].getX(), right[2].getY());
            dr = std::min(d3, std::min(d1, d2));
            if (dr == d2) {
                return Distance(dr, right[1], right[2]);
            }
            if (dr == d3) {
                return Distance(dr, right[0], right[2]);
            }
        }
        dr = d1;
        return Distance(dr, right[0], right[1]);
    }

    // Find the middle point
    int mid = v.size() / 2;
    Point midPoint = v[mid];

    // Consider the vertical line passing through the middle point calculate
    // the smallest distance dl on left of middle point and dr on right side
    Distance dl = recursive(left);
    Distance dr = recursive(right);

    // Find the smaller of two distances
    Distance d = dl.compare(dl, dr);

    // Build an array strip[] that contains points close (closer than d)
    // to the line passing through the middle point
    std::vector<Point> strip;
    for (size_t i = 0; i < v.size(); i++)
        if (abs(v[i].getX() - midPoint.getX()) < d.getdist())
            strip.push_back(v[i]);

    strip = Point::getSortedY(strip);
    double temp;

    for (size_t p1 = 0; p1 != strip.size(); p1++)
        for (size_t p2 = p1 + 1; p2 != p1 + 15; p2++) {
            temp = distance(strip[p1].getX(), strip[p1].getY(), strip[p2].getX(), strip[p2].getY());
            if (temp < d.getdist())
                d = Distance(temp, strip[p1], strip[p2]);
        }
    // Find the closest points in strip. Return the minimum of d and closest distance is strip[]
    return d;
}

Distance part3() {
    std::vector<Point> v = {};

    //pulling the numbers from the file
    std::vector<std::string> vec = {};
    std::string line;
    std::ifstream points("points.txt");

    if (points.is_open()) {
        while (getline(points, line)) {
            std::stringstream ss(line);
            while (getline(ss, line, ' ')) {
                if (line.length() > 2)
                    vec.push_back(line);
            }
        }
        points.close();
    }
    else
        std::cout << "Unable to open file";

    std::vector<std::string>::iterator it;
    for (it = vec.begin(); it != vec.end(); it++) {
        std::string s1 = *it;
        it++;
        std::string s2 = *it;
        v.push_back(Point(std::stod(s1), std::stod(s2)));
    }

    std::clock_t start = std::clock();

    v = Point::getSorted(v);

    Distance d = recursive3(v);
    double duration = ((double)std::clock() - (double)start) / (double)CLOCKS_PER_SEC;
    duration3 = duration;

    return d;
}

void knuth(std::vector<Point> &v) {
    std::vector<Point> knuth;
    srand(time(NULL));
    int index = 0;
    while (index < v.size()) {
        double random = double(rand()) / RAND_MAX;
        int newindex = (int)(index + (random * (v.size() - index)));
        Point temp = v[index];
        v[index] = v[newindex];
        v[newindex] = temp;
        index++;
        //index = rand() % v.size();

        /*std::cout << index << " ";
        v[index].print();
        std::cout << "\n";*/
        //std::cout << v.size() << "\n";
        
        //knuth.push_back(v[index]);
        //v.erase(v.begin() + index);
    }
    //return knuth;
}

Distance returndist(Point a, Point b) {
    return Distance(distance(a.getX(), a.getY(), b.getX(), b.getY()), a, b);
}

/*std::unordered_map<Box, Point, MyHashFunction> makedict(std::vector<Point> v, double sidelength) {
    std::unordered_map<Box, Point, MyHashFunction> map;
    std::cout << "@4\n";
    for (double i = 0; i <= 1; i += sidelength) {
        for (double j = 0; j <= 1; j += sidelength) {
            for (size_t index = 0; index < v.size(); index++) {
                if (v[index].getX() < i + sidelength && v[index].getY() < j + sidelength) {
                    //std::cout << i / sidelength << " " << j / sidelength << "\t" << v[index].getX() << " " << v[index].getY() << "\n";
                    Box b(i / sidelength, j / sidelength);
                    map[b] = v[index];
                    v.erase(v.begin() + index);
                }
            }
        }
    }
    return map;
}*/

Box pointtobox(Point p, Distance smallestdist) {
    unsigned long long int xcoord = (unsigned long long int)(p.getX() / smallestdist.getdist() / 2);
    unsigned long long int ycoord = (unsigned long long int)(p.getY() / smallestdist.getdist() / 2);
    //Box(xcoord, ycoord).print();
    //std::cout << "\n";
    return Box(xcoord, ycoord);
}

std::unordered_map<Box, Point, MyHashFunction> makedict(std::unordered_map<Box, Point, MyHashFunction> map, Distance smallestdist) {
    std::unordered_map<Box, Point, MyHashFunction> newmap;
    for (auto p : map)
        newmap[pointtobox(p.second, smallestdist)] = p.second;
    /*//for (double i = 0; i <= 1; i += sidelength) {
    //    for (double j = 0; j <= 1; j += sidelength) {
            for (size_t index = 0; index < v.size(); index++) {
                if (v[index].getX() < i + sidelength && v[index].getY() < j + sidelength) {
                    //std::cout << i / sidelength << " " << j / sidelength << "\t" << v[index].getX() << " " << v[index].getY() << "\n";
                    Box b(i / sidelength, j / sidelength);
                    map[b] = v[index];
                    v.erase(v.begin() + index);
                }
            }
    //    }
    //}*/
    return newmap;
}

Distance algorithm(std::vector<Point> &v) {
    Distance smallestdist = returndist(v[0], v[1]);
    std::unordered_map<Box, Point, MyHashFunction> map;
    map[pointtobox(v[0], smallestdist)] = v[0];
    map[pointtobox(v[1], smallestdist)] = v[1];

    for (size_t i = 2; i < v.size(); i++) {
        //Box b(v[i].getX() / smallestdist.getdist() / 2, v[i].getY() / smallestdist.getdist() / 2);
        Box b(pointtobox(v[i], smallestdist));
        //std::cout << b.getX() << " " << b.getY() << " ";

        //std::unordered_map<Box, Point, MyHashFunction>::iterator it;

        Distance mindist25 = smallestdist;
        Point p1, p2;

        /*for (auto const& it : map) {
            if (it.first == b) {
                smallestdist = returndist(it.second, v[i]);
                std::unordered_map<Box, Point, MyHashFunction> temp = makedict(map, smallestdist);
                map = temp;
                temp.clear();
                break;
            }
        }*/
        //map[pointtobox(v[i], smallestdist)] = v[i];

        unsigned long long int x = b.getX();
        unsigned long long int y = b.getY();
        unsigned long long int minx, maxx, miny, maxy;
        minx = x - 2;
        maxx = x + 2;
        miny = y - 2;
        maxy = y + 2;

        /*unsigned long long int numboxes = 1 / smallestdist.getdist() / 2;
        if (x < 2)
            minx = 0;
        else
            minx = x - 2;

        if (x + 2 > numboxes - 1)
            maxx = numboxes - 1;
        else
            maxx = x + 2;

        if (y - 2 < 0)
            miny = 0;
        else
            miny = y - 2;

        if (y + 2 > numboxes - 1)
            maxy = numboxes - 1;
        else
            maxy = y + 2;*/

        //std::vector<Distance> smallerdist = {};
        for (unsigned long long int i = minx; i <= maxx; i++) {
            for (unsigned long long int j = miny; j <= maxy; j++) {
                //if (i != x && j != y) {
                if (returndist(v[i], map[Box(i, j)]).getdist() < mindist25.getdist())
                    mindist25 = returndist(v[i], map[Box(i, j)]);
                //}
            }
        }
        if (mindist25.getdist() >= smallestdist.getdist())
            map[b] = v[i];
        else {
            // redo the map
            std::unordered_map<Box, Point, MyHashFunction> newmap;
            for (auto p : map)
                newmap[pointtobox(p.second, smallestdist)] = p.second;
            //std::unordered_map<Box, Point, MyHashFunction> temp = makedict(map, mindist25);
            map = newmap;
            newmap.clear();
        }

        /*for (auto a : smallerdist) {
            if (a.getdist() < smallestdist.getdist())
                smallestdist = a;
        }*/
        /*std::unordered_map<Box, Point, MyHashFunction> temp = makedict(map, smallestdist);
        map = temp;
        temp.clear();*/
    }

    return smallestdist;
    /*Distance S = returndist(v[0], v[1]);
    Distance d;

    std::unordered_map<Box, Point, MyHashFunction> map;
    map = makedict(v, S.getdist() / 2);
    for (size_t i = 0; i < v.size(); i++) {
        // detrmine the subsquare containing pi
        Box subsquare;
        for (auto e : map)
            if (e.second.getX() == v[i].getX() && e.second.getY() == v[i].getY())
                subsquare = e.first;

        std::cout << v[i].getX() << " " << v[i].getY() << "\n";
        std::cout << subsquare.getX() << " " << subsquare.getY() << "\n";

        // look up the 25 subsquares close to pi
        unsigned long long int lowerleftX;
        unsigned long long int lowerleftY;
        unsigned long long int upperrightX;
        unsigned long long int upperrightY;
        if (subsquare.getX() >= 2)
            lowerleftX = subsquare.getX() - 2;
        else
            lowerleftX = 0;
        if (subsquare.getY() >= 2)
            lowerleftY = subsquare.getY() - 2;
        else
            lowerleftY = 0;
        if (subsquare.getX() < 1 / S.getdist() - 2)
            upperrightX = subsquare.getX() - 2;
        else
            upperrightX = 1 / S.getdist();
        if (subsquare.getY() < 1 / S.getdist() - 2)
            upperrightY = subsquare.getY() - 2;
        else
            upperrightY = 0;

        Box start(lowerleftX, lowerleftY);
        Box end(upperrightX, upperrightY);

        bool inbox = false;

        // Computer the distance from pi to any points founds in these subsquares
        for (size_t j = 0; j < v.size(); j++) {
            for (unsigned long long int i = lowerleftX; i <= upperrightX; i++) {
                for (unsigned long long int j = lowerleftY; j <= upperrightY; j++) {
                    Box b(i, j);
                    if (map[b].same(v[j])) {
                        inbox = true;
                        //std::cout << "@6\n";
                    }
                }
            }
            if (inbox == true) {
                d = returndist(v[i], v[j]);
                if (d.getdist() < S.getdist()) {
                    map.clear();
                    map = makedict(v, d.getdist() / 2);
                    for (size_t k = 0; k <= i; k++) {
                        for (auto e : map)
                            if (e.second.getX() == v[k].getX() && e.second.getY() == v[k].getY())
                                subsquare = e.first;
                    }
                }
                else {
                    // put in pi into dictionary
                    map[subsquare] = v[i];
                }
            }
            inbox = false;
        }
    }
    return d;*/
}

Distance part4() {
    std::vector<Point> v = {};

    //pulling the numbers from the file
    std::vector<std::string> vec = {};
    std::string line;
    std::ifstream points("points.txt");

    if (points.is_open()) {
        while (getline(points, line)) {
            std::stringstream ss(line);
            while (getline(ss, line, ' ')) {
                if (line.length() > 2)
                    vec.push_back(line);
            }
        }
        points.close();
    }
    else
        std::cout << "Unable to open file";

    std::vector<std::string>::iterator it;
    for (it = vec.begin(); it != vec.end(); it++) {
        std::string s1 = *it;
        it++;
        std::string s2 = *it;
        v.push_back(Point(std::stod(s1), std::stod(s2)));
    }

    std::clock_t start = std::clock();

    // knuth shuffle
    knuth(v);

    Distance d = algorithm(v);
    double duration = ((double)std::clock() - (double)start) / (double)CLOCKS_PER_SEC;
    duration4 = duration;

    return d;
}

int main() {
    results.open("results.txt");

    //Distance d1 = part1();
    //results << std::fixed;
    //results << std::setprecision(23) << d1.getdist() << "\t";
    //results << "(" << d1.getp1().getX() << " , " << d1.getp2().getY() << ")";
    //std::cout << "  ";
    //results << "(" << d1.getp1().getX() << " , " << d1.getp2().getY() << ")";
    //results << "\nduration: " << duration1 << " seconds";
    // d1.print();
    //std::cout << "\nduration: " << duration1 << " seconds\n";

    //Distance d2 = part2();
    //results << std::fixed;
    //results << std::setprecision(23) << d2.getdist() << "\t";
    //results << "(" << d2.getp1().getX() << " , " << d2.getp1().getY() << ")";
    //std::cout << "  ";
    //results << "(" << d2.getp2().getX() << " , " << d2.getp2().getY() << ")";
    //results << "\nduration: " << duration2 << " seconds";
    //d2.print();
    //std::cout << "\nduration: " << duration2 << " seconds\n";

    Distance d3 = part3();
    results << std::fixed;
    results << std::setprecision(23) << d3.getdist() << "\t";
    results << "(" << d3.getp1().getX() << " , " << d3.getp1().getY() << ") ";
    //std::cout << "  ";
    results << "(" << d3.getp2().getX() << " , " << d3.getp2().getY() << ")";
    results << "\nduration: " << duration3 << " seconds";

    d3.print();
    std::cout << "\nduration: " << duration3 << " seconds\n";

    Distance d4 = part4();
    results << std::fixed;
    results << std::setprecision(23) << "\n" << d4.getdist() << "\t";
    results << "(" << d4.getp1().getX() << " , " << d4.getp1().getY() << ")";
    std::cout << "  ";
    results << "(" << d4.getp2().getX() << " , " << d4.getp2().getY() << ")";
    results << "\nduration: " << duration4 << " seconds";

    d4.print();
    std::cout << "\nduration: " << duration4 << " seconds\n";

    results.close();
}