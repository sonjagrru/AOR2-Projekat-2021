/**
 * @file vigenere_cipher.cpp
 * @brief Implementation of [Vigenère cipher](https://en.wikipedia.org/wiki/Vigen%C3%A8re_cipher) algorithm.
 *
 * @details
 * The Vigenère cipher is a method of encrypting alphabetic text by using a series of interwoven vigenere
 * ciphers, based on the letters of a keyword. It employs a form of polyalphabetic substitution.
 *
 * ### Algorithm
 * The encryption can also be represented using modular arithmetic by first transforming
 * the letters into numbers, according to the scheme, A → 0, B → 1, ..., Z → 25.
 * Encryption of \f$i^{th}\f$ character in Message M by key K can be described mathematically as,
 *
 * \f[ E_{K}(M_{i}) = (M_{i} + K_{i})\;\mbox{mod}\; 26\f]
 *
 * while decryption of \f$i^{th}\f$ character in Cipher C by key K can be described mathematically as,
 *
 * \f[ D_{k}(C_{i}) = (C_{i} - K_{i} + 26)\;\mbox{mod}\; 26\f]
 *
 * Where \f$K_{i}\f$ denotes corresponding character in key. If \f$|key| < |text|\f$ than
 * same key is repeated untill their lengths are equal.
 *
 * For Example,
 * If M = "ATTACKATDAWN" and K = "LEMON" than K becomes "LEMONLEMONLE".
 *
 * \note Rather than creating new key of equal length this program does this by using modular index for key
 * (i.e. \f$(j + 1) \;\mbox{mod}\; |\mbox{key}|\f$)
 *
 * \note This program implements Vigenère cipher for only uppercase English alphabet characters (i.e. A-Z).
 *
 * @author [Deep Raval](https://github.com/imdeep2905)
 */
#include <iostream>
#include <string>
#include <cassert>
#include <immintrin.h>

#include "emmintrin.h"
#include "_Timer.h"

using namespace std;

 /** \namespace ciphers
  * \brief Algorithms for encryption and decryption
  */
namespace ciphers {
    /** \namespace vigenere
     * \brief Functions for [vigenère cipher](https://en.wikipedia.org/wiki/Vigen%C3%A8re_cipher) algorithm.
     */
    namespace vigenere {
        namespace {
            /**
             * This function finds character for given value (i.e.A-Z)
             * @param x value for which we want character
             * @return  corresponding character for perticular value
             */
            inline char get_char(const int x) {
                // By adding 65 we are scaling 0-25 to 65-90. 
                // Which are in fact ASCII values of A-Z. 
                return char(x + 65);
            }
            /**
             * This function finds value for given character (i.e.0-25)
             * @param c character for which we want value
             * @return returns corresponding value for perticular character
             */
            inline int get_value(const char c) {
                // A-Z have ASCII values in range 65-90.
                // Hence subtracting 65 will scale them to 0-25.
                return int(c - 65);
            }
        } // Unnamed namespace
        /**
         * Encrypt given text using vigenere cipher.
         * @param text text to be encrypted
         * @param key to be used for encryption
         * @return new encrypted text
         */
        std::string encrypt(const std::string& text, const std::string& key) {
            std::string encrypted_text = ""; // Empty string to store encrypted text
            // Going through each character of text and key
            // Note that key is visited in circular way hence  j = (j + 1) % |key|
            for (size_t i = 0, j = 0; i < text.length(); i++, j = (j + 1) % key.length()) {
                int place_value_text = get_value(text[i]); // Getting value of character in text
                int place_value_key = get_value(key[j]); // Getting value of character in key
                place_value_text = (place_value_text + place_value_key) % 26; // Applying encryption
                char encrypted_char = get_char(place_value_text); // Getting new character from encrypted value
                encrypted_text += encrypted_char; // Appending encrypted character
            }
            return encrypted_text; // Returning encrypted text
        }

        /*OPTIMIZOVANO*/
        std::string encryptO(const std::string& text, const std::string& key) {
            std::string encrypted_text = ""; 

            _mm_prefetch((char*)(&text[0]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[1]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[2]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[3]), _MM_HINT_T1);

            size_t i = 4, j = 0;
            for (; i < text.length(); i++, j = (j + 1) % key.length()) {
                _mm_prefetch((char*)(&text[i]), _MM_HINT_T1);
                int place_value_text = get_value(text[i-4]); 
                int place_value_key = get_value(key[j]); 
                place_value_text = (place_value_text + place_value_key) % 26; 
                char encrypted_char = get_char(place_value_text); 
                encrypted_text += encrypted_char; 
            }
            for (i = 0; i < 4; i++)
            {
                int place_value_text = get_value(text[text.length() - 4 + i]);
                int place_value_key = get_value(key[j]);
                place_value_text = (place_value_text + place_value_key) % 26;
                char encrypted_char = get_char(place_value_text);
                encrypted_text += encrypted_char;
                j = (j + 1) % key.length();
            }

            return encrypted_text; 
        }
        /**
         * Decrypt given text using vigenere cipher.
         * @param text text to be decrypted
         * @param key key to be used for decryption
         * @return new decrypted text
         */
        std::string decrypt(const std::string& text, const std::string& key) {
            // Going through each character of text and key
            // Note that key is visited in circular way hence  j = (j + 1) % |key|
            std::string decrypted_text = ""; // Empty string to store decrypted text
            for (size_t i = 0, j = 0; i < text.length(); i++, j = (j + 1) % key.length()) {
                int place_value_text = get_value(text[i]); // Getting value of character in text
                int place_value_key = get_value(key[j]); // Getting value of character in key
                place_value_text = (place_value_text - place_value_key + 26) % 26; // Applying decryption
                char decrypted_char = get_char(place_value_text); // Getting new character from decrypted value
                decrypted_text += decrypted_char; // Appending decrypted character
            }
            return decrypted_text; // Returning decrypted text
        }

        /*OPTIMIZOVANO PREFETCH*/
        std::string decryptO(const std::string& text, const std::string& key) {

            _mm_prefetch((char*)(&text[0]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[1]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[2]), _MM_HINT_T1);
            _mm_prefetch((char*)(&text[3]), _MM_HINT_T1);

            std::string decrypted_text = ""; // Empty string to store decrypted text

            size_t i = 4, j = 0;
            for (; i < text.length(); i++, j = (j + 1) % key.length()) {
                _mm_prefetch((char*)(&text[i]), _MM_HINT_T1);
                int place_value_text = get_value(text[i-4]);
                int place_value_key = get_value(key[j]); 
                place_value_text = (place_value_text - place_value_key + 26) % 26; 
                char decrypted_char = get_char(place_value_text); 
                decrypted_text += decrypted_char; 
            }
            for (i = 0; i < 4; i++)
            {
                int place_value_text = get_value(text[text.length() - 4 + i]);
                int place_value_key = get_value(key[j]);
                place_value_text = (place_value_text - place_value_key + 26) % 26;
                char decrypted_char = get_char(place_value_text);
                decrypted_text += decrypted_char;
                j = (j + 1) % key.length();
            }
            
            return decrypted_text; 
        }
    } // namespace vigenere
} // namespace ciphers

/**
 * Function to test above algorithm
 */
void test() {
    // Test 1
    std::string text1 = "NIKOLATESLA";
    std::string encrypted1, decrypted1;
    StartTimer(ORIGINAL)
    encrypted1 = ciphers::vigenere::encrypt(text1, "TESLA");
    decrypted1 = ciphers::vigenere::decrypt(encrypted1, "TESLA");
    EndTimer
    assert(text1 == decrypted1);
    std::cout << "Original text : " << text1;
    std::cout << " , Encrypted text (with key = TESLA) : " << encrypted1;
    std::cout << " , Decrypted text : " << decrypted1 << std::endl;
    
    
    //OPTIMIZOVANOPREFETCH
    text1 = "NIKOLATESLA";
    StartTimer(OPTIMIZOVANOPREFETCH)
    encrypted1 = ciphers::vigenere::encryptO(text1, "TESLA");
    decrypted1 = ciphers::vigenere::decryptO(encrypted1, "TESLA");
    EndTimer
    //assert(text1 == decrypted1);
    std::cout << "Original text : " << text1;
    std::cout << " , Encrypted text (with key = TESLA) : " << encrypted1;
    std::cout << " , Decrypted text : " << decrypted1 << std::endl;

    // Test 2
    std::string text2 = "";
    for (int i = 0; i < 10000; i++)
    {
        int a = rand()%26;
        text2 += ciphers::vigenere::get_char(a);
    }
    std::string encrypted2, decrypted2;
    StartTimer(ORIGINAL)
    encrypted2 = ciphers::vigenere::encrypt(text2, "REALLY");
    decrypted2 = ciphers::vigenere::decrypt(encrypted2, "REALLY");
    EndTimer
    //assert(text2 == decrypted2);
    std::cout << "Original text : " << text2;
    std::cout << " ,\n Encrypted text (with key = REALLY) : " << encrypted2;
    std::cout << " ,\n Decrypted text : " << decrypted2 << std::endl;

    StartTimer(OPTIMIZOVANOPREFETCH)
    encrypted2 = ciphers::vigenere::encrypt(text2, "REALLY");
    decrypted2 = ciphers::vigenere::decrypt(encrypted2, "REALLY");
    EndTimer
    //assert(text2 == decrypted2);
    std::cout << "Original text : " << text2;
    std::cout << " ,\n Encrypted text (with key = REALLY) : " << encrypted2;
    std::cout << " ,\n Decrypted text : " << decrypted2 << std::endl;
}

/** Driver Code */
int main() {
    // Testing
    test();
    return 0;
}
