#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

struct Pair
{
    explicit Pair() {};
    explicit Pair(int i, int j) : f_coulomb(3), f_vw(3), f_bond(3), f_angle(3), f_torsion(3), fij(3) {
            atomi_id = i;
            atomj_id = j;
    };

    explicit Pair(int i, int j, std::vector<double> f_coulomb, std::vector<double> f_vw, std::vector<double> f_bond, std::vector<double> f_angle, std::vector<double> f_torsion) : \
        f_coulomb(3), f_vw(3), f_bond(3), f_angle(3), f_torsion(3), fij(3){
        this->f_coulomb = f_coulomb;
        this->f_vw = f_vw;
        this->f_bond = f_bond;
        this->f_angle = f_angle;
        this->f_torsion = f_torsion;
    };

    std::vector<double> &get_fij() {
        std::vector<double> force(3);
        for (size_t i = 0; i < this->f_coulomb.size(); i++)
        {
            force[i] = f_coulomb[i] + f_vw[i] + f_bond[i] + \
                f_angle[i] + f_torsion[i];
        }
        this->fij = force;
        return fij;
    };

private:
    int atomi_id;
    int atomj_id;
    std::vector<double> f_coulomb;
    std::vector<double> f_vw;
    std::vector<double> f_bond;
    std::vector<double> f_angle;
    std::vector<double> f_torsion;
    std::vector<double> fij;
};

struct Pairs
{
    explicit Pairs(int num_pairs) : _pairs(num_pairs) {
        this->_iter_id = 0;
        this->_num_pairs = num_pairs;
    };

    void add_pair(int atom1_id, int atom2_id, std::vector<double> f_coulomb, std::vector<double> f_vw, std::vector<double> f_bond, \
                std::vector<double> f_angle, std::vector<double> f_torsion) {
        Pair pair(atom1_id, atom2_id, f_coulomb, f_angle, f_bond, f_angle, f_torsion);
        this->_pairs[this->_iter_id];
        this->_iter_id ++;
    };

    void add_pairs(std::vector<Pair> pairs) {
        this->_pairs = pairs;
    };

    void update_coulombic_force(int atom1_id, int atom2_id, std::vector<double> fij) {
    };

    void reset() {
        this->_iter_id = 0;
    };

    int num_pairs() {
        return this->_num_pairs;
    }

private:
    std::vector<Pair> _pairs;
    int _iter_id;
    int _num_pairs;
};