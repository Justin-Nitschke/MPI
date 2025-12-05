struct Node {
    float com_x, com_y;
    float total_mass;
    Node *children[4];
    float square_x, square_y;
    float half_width;
};