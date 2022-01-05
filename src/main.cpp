#include <mprocessor/music.h>

void basscallback()
{
    cout << "bass detected" << endl;
}

/**
 * @brief new comment
 * 
 * @return int 
 */
int main()
{
    init(&basscallback);  
    while(1)
        update();

    return 0;
}