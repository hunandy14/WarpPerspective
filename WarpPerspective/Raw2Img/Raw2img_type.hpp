/*****************************************************************
Name :
Date : 2017/06/23
By   : CharlotteHonG
Final: 2017/06/23
*****************************************************************/
#pragma once
/*
    ########   ######   ########
    ##     ## ##    ##  ##     ##
    ##     ## ##        ##     ##
    ########  ##   #### ########
    ##   ##   ##    ##  ##     ##
    ##    ##  ##    ##  ##     ##
    ##     ##  ######   ########
*/

enum RGBs {R, G, B};
class RGBs_t {
public:
    RGBs_t(size_t rgb): rgb(RGBs(rgb)) {
        if(rgb>2) {throw std::range_error("range only 0~2");}
    }
    inline operator RGBs() {return rgb;}
private:
    RGBs rgb;
};
/*
    ######## #### ##       ########         ##     ##
    ##        ##  ##       ##               ##     ##
    ##        ##  ##       ##               ##     ##
    ######    ##  ##       ######           #########
    ##        ##  ##       ##               ##     ##
    ##        ##  ##       ##               ##     ##
    ##       #### ######## ######## ####### ##     ##
*/

// 檔案檔頭 (BITMAPFILEHEADER)
#pragma pack(2) // 調整對齊
struct BmpFileHeader{
    uint16_t bfTybe=0x4d42;
    uint32_t bfSize;
    uint16_t bfReserved1=0;
    uint16_t bfReserved2=0;
    uint32_t bfOffBits=54;
    void pri(){
        using namespace std;
        cout << "BmpFileHeader" << endl;
        cout << "  bfTybe=" << bfTybe << endl;
        cout << "  bfSize=" << bfSize << endl;
        cout << "  bfReserved1=" << bfReserved1 << endl;
        cout << "  bfReserved2=" << bfReserved2 << endl;
        cout << "  bfOffBits=" << bfOffBits << endl;
    }
    // fstream
    friend std::ofstream& operator>>(
        std::ofstream& is, BmpFileHeader& obj);
    friend std::ifstream& operator<<(
        std::ifstream& os, const BmpFileHeader& obj);
    // ostream
    friend std::ostream& operator<<(
        std::ostream& os, const BmpFileHeader& obj);
};
#pragma pack() // 恢復對齊為預設
inline std::ofstream& operator<<(
    std::ofstream& os, const BmpFileHeader& obj){
    os.write((char*)&obj, sizeof(obj));
    return os;
}
inline std::ifstream& operator>>(
    std::ifstream& is, BmpFileHeader& obj){
    is.read((char*)&obj, sizeof(obj));
    return is;
}
inline std::ostream& operator<<(
    std::ostream& os, const BmpFileHeader& obj){
    using std::cout;
    using std::endl;
    cout << "# BmpFileHeader" << endl;
    cout << "    bfTybe      = " << obj.bfTybe << endl;
    cout << "    bfSize      = " << obj.bfSize << endl;
    cout << "    bfReserved1 = " << obj.bfReserved1 << endl;
    cout << "    bfReserved2 = " << obj.bfReserved2 << endl;
    cout << "    bfOffBits   = " << obj.bfOffBits;
    return os;
}
/*
    #### ##    ## ########  #######          ##     ##
     ##  ###   ## ##       ##     ##         ##     ##
     ##  ####  ## ##       ##     ##         ##     ##
     ##  ## ## ## ######   ##     ##         #########
     ##  ##  #### ##       ##     ##         ##     ##
     ##  ##   ### ##       ##     ##         ##     ##
    #### ##    ## ##        #######  ####### ##     ##
*/

// 圖片資訊 (BITMAPINFOHEADER)
#pragma pack(2) // 調整對齊
struct BmpInfoHeader{
    uint32_t biSize=40;
    uint32_t biWidth;
    uint32_t biHeight;
    uint16_t biPlanes=1; // 1=defeaul, 0=custom
    uint16_t biBitCount;
    uint32_t biCompression=0;
    uint32_t biSizeImage;
    uint32_t biXPelsPerMeter=0; // 72dpi=2835, 96dpi=3780
    uint32_t biYPelsPerMeter=0; // 120dpi=4724, 300dpi=11811
    uint32_t biClrUsed=0;
    uint32_t biClrImportant=0;
    // fstream
    friend std::ofstream& operator>>(
        std::ofstream& is, BmpInfoHeader& obj);
    friend std::ifstream& operator<<(
        std::ifstream& os, const BmpInfoHeader& obj);
    // ostream
    friend std::ostream& operator<<(
        std::ostream& os, const BmpInfoHeader& obj);
};
#pragma pack() // 恢復對齊為預設
inline std::ofstream& operator<<(
    std::ofstream& os, const BmpInfoHeader& obj){
    os.write((char*)&obj, sizeof(obj));
    return os;
}
inline std::ifstream& operator>>(
    std::ifstream& is, BmpInfoHeader& obj){
    is.read((char*)&obj, sizeof(obj));
    return is;
}
inline std::ostream& operator<<(
    std::ostream& os, const BmpInfoHeader& obj){
    using std::cout;
    using std::endl;
    cout << "# BmpInfoHeader" << endl;
    cout << "    biSize          = " << obj.biSize << endl;
    cout << "    biWidth         = " << obj.biWidth << endl;
    cout << "    biHeight        = " << obj.biHeight << endl;
    cout << "    biPlanes        = " << obj.biPlanes << endl;
    cout << "    biBitCount      = " << obj.biBitCount << endl;
    cout << "    biCompression   = " << obj.biCompression << endl;
    cout << "    biSizeImage     = " << obj.biSizeImage << endl;
    cout << "    biXPelsPerMeter = " << obj.biXPelsPerMeter << endl;
    cout << "    biYPelsPerMeter = " << obj.biYPelsPerMeter << endl;
    cout << "    biClrUsed       = " << obj.biClrUsed << endl;
    cout << "    biClrImportant  = " << obj.biClrImportant;
    return os;
}
