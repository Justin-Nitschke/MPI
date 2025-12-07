#include "Body.hpp"

struct Node {
    float com_x, com_y;
    float total_mass;
    Body *body;
    Node *children[4];
    float square_x, square_y;
    float half_width;
    bool is_internal;

    Node(float x, float y, float hw) : square_x(x), square_y(y), half_width(hw) {
        body = nullptr;
        is_internal = false;
        total_mass = 0;
        com_x = 0;
        com_y = 0;
        for (int i = 0; i < 4; i++) children[i] = nullptr;
    }

    ~Node() {
        for (int i = 0; i < 4; i++) {
            if (children[i]) delete children[i];
        }
    }

    void insert(Body *b) {
        if (is_internal) {
            children[get_quadrant(b)]->insert(b);
        } else if (!body) {
            body = b;
        } else {
            if (body->x_position == b->x_position && body->y_position == b->y_position) {
                body->mass += b->mass;
                return;
            }

            Body *old_b = body;

            is_internal = true;
            body = nullptr;

            split();

            children[get_quadrant(old_b)]->insert(old_b);
            children[get_quadrant(b)]->insert(b);
        }
    }

    void center_of_mass() {
        if (!body && !is_internal) return;

        if (!is_internal) {
            total_mass = body->mass;
            com_x = body->x_position;
            com_y = body->y_position;
            return;
        }

        float mass_sum = 0;
        float moment_x = 0;
        float moment_y = 0;

        for (int i = 0; i < 4; i++) {
            if (children[i]) {
                children[i]->center_of_mass();

                float m = children[i]->total_mass;
                if (m > 0) {
                    mass_sum += m;
                    moment_x += m * children[i]->com_x;
                    moment_y += m * children[i]->com_y;
                }
            }
        }

        total_mass = mass_sum;
        if(total_mass > 0) {
            com_x = moment_x / total_mass;
            com_y = moment_y / total_mass;
        } else {
            com_x = 0;
            com_y = 0;
        }
    }

    void force(Body *b, double theta) {
        if (b == body || total_mass == 0) {
            return;
        }

        float dx = com_x - b->x_position;
        float dy = com_y - b->y_position;
        float dist = std::sqrt(dx * dx + dy * dy);
        
        float rlimit = 0.03;
        dist = std::max(rlimit, dist);

        float width = 2 * half_width;

        if (!is_internal || (width / dist) < theta) {
            float G = 0.1;
            float force = (G * b->mass * total_mass) / (dist * dist);
            
            float fx = force * (dx / dist);
            float fy = force * (dy / dist);

            b->x_force += fx;
            b->y_force += fy;
        } else {
            for (int i = 0; i < 4; i++) {
                if (children[i]) {
                    children[i]->force(b, theta);
                }
            }
        }
    }

    void split() {
        float quarter = half_width / 2;
        float center_x = square_x + half_width;
        float center_y = square_y + half_width;

        children[0] = new Node(square_x, center_y, quarter);
        children[1] = new Node(center_x, center_y, quarter);
        children[2] = new Node(square_x, square_y, quarter);
        children[3] = new Node(center_x, square_y, quarter);
    }

    int get_quadrant(Body *b) {
        double center_x = square_x + half_width;
        double center_y = square_y + half_width;

        bool left = b->x_position < center_x;
        bool bottom = b->y_position < center_y;

        if (left && !bottom) return 0;
        else if (!left && !bottom) return 1;
        else if (left && bottom) return 2;
        else return 3;
    }
};