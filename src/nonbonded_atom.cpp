#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

class Atom
{
public:
    void add_meta_properties(std::string name, std::string element, int id, int residue_id)
    {
        this->name = name;
        this->element = element;
        this->id = id;
        this->residue_id = residue_id;
    };

    void add_nonbonded_params(double mass, double charge, double sigma, double epsilon)
    {
        this->mass = mass;
        this->charge = charge;
        this->sigma = sigma;
        this->epsilon = epsilon;
    };

    void update_position(std::vector<double> position)
    {
        this->position = position;
    };

    void get_params() const
    {
        std::cout << "------------------------------" << std::endl;
        std::cout << "name     : " << name << std::endl;
        std::cout << "element  : " << element << std::endl;
        std::cout << "id       : " << id << std::endl;
        std::cout << "residue_id  : " << residue_id << std::endl;
        std::cout << "charge   : " << charge << std::endl;
        std::cout << "sigma    : " << sigma << std::endl;
        std::cout << "epsilon  : " << epsilon << std::endl;
        std::cout << "position : (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
        std::cout << "------------------------------" << std::endl
                  << std::endl;
    };

private:
    std::string name;
    std::string element;
    int id;
    int residue_id;
    double mass;
    double charge;
    double sigma;
    double epsilon;
    std::vector<double> position;

public:
    std::string get_name()
    {
        return this->name;
    };

    std::string get_element()
    {
        return this->element;
    };

    int get_id()
    {
        return this->id;
    };

    int get_residue_id()
    {
        return this->residue_id;
    };

    double get_charge()
    {
        return this->charge;
    };

    double get_sigma()
    {
        return this->sigma;
    };

    double get_epsilon()
    {
        return this->epsilon;
    };

    std::vector<double> get_position()
    {
        return this->position;
    };
};