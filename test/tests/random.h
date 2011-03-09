dmatrix m1(4, 4);
dvector v1(4);

juzhen::UpdateRandomSeed();

std::cout << "juzhen::RandVector<double>(5, 1) = " << juzhen::RandVector<double>(5, 1) << std::endl;
std::cout << "juzhen::RandVector<double>(5) = " << juzhen::RandVector<double>(5) << std::endl;
std::cout << "juzhen::RandMatrix<double>(5, 5, 1) = " << juzhen::RandMatrix<double>(5, 5, 1) << std::endl;
std::cout << "juzhen::RandMatrix<double>(5, 5) = " <<juzhen::RandMatrix<double>(5, 5) << std::endl;

std::cout << "Randomize(m1, 10) = \n" << Randomize(m1, 10) << std::endl;
std::cout << "Randomize(m1) = \n" << Randomize(m1) << std::endl;
std::cout << "Randomize(v1, 10) = \n" << Randomize(v1, 10) << std::endl;
std::cout << "Randomize(v1) = \n" << Randomize(v1) << std::endl;

//std::cout << RandMatrix(5, 5, 1) << std::endl;
