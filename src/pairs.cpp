#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

class Pair {
public:
    Pair(int i, int j) : coulombic_force(3), vander_walls_force(3), bonded_force(3), \
        angle_force_contrib(3), torsion_force_contrib(3), fij(3){
            atomi_id = i;
            atomj_id = j;
    };

    std::vector<double> &get_fij() {
        std::vector<double> force(3);
        for (size_t i = 0; i < this->coulombic_force.size(); i++)
        {
            force[i] = coulombic_force[i] + vander_walls_force[i] + bonded_force[i] + \
                angle_force_contrib[i] + torsion_force_contrib[i];
        }
        this->fij = force;
        return fij;
    };

private:
    int atomi_id;
    int atomj_id;
    std::vector<double> coulombic_force;
    std::vector<double> vander_walls_force;
    std::vector<double> bonded_force;
    std::vector<double> angle_force_contrib;
    std::vector<double> torsion_force_contrib;
    std::vector<double> fij;
};

class Pairs
{
public:
    void add_pair(int atom1_id, int atom2_id) {
        Pair pair(atom1_id, atom2_id);
        this->pairs.push_back(pair);
    };

    void add_pairs(std::vector<Pair> pairs) {
        this->pairs = pairs;
    };

    void update_coulombic_force(int atom1_id, int atom2_id, std::vector<double> fij) {
    };

private:
    std::vector<Pair> pairs;
};