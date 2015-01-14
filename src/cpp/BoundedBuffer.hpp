// Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following
// disclaimer in the documentation and/or other materials provided
// with the distribution.
//
// * Neither the name of Pacific Biosciences nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.



///
/// Templated, thread-safe buffer container, uses uses boost::circular buffer 
/// bounded by a given capacity specified by the caller.  When the buffer is 
/// full, the push waits for an open spot.  When the buffer is empty, the pop
/// waits for an item to be present.  Condition variables are used to signal
/// the state of the buffer.
///
template <class T>
class BoundedBuffer {
public:
    typedef boost::circular_buffer<T> container_type;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::value_type value_type;
    typedef typename boost::call_traits<value_type>::param_type param_type;

    explicit BoundedBuffer(size_type capacity) : unread_(0), container_(capacity) {}

    void push(param_type item) {
        boost::mutex::scoped_lock lock(mutex_);
        not_full_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_full, this));
        container_.push_front(item);
        ++unread_;
        lock.unlock();
        not_empty_.notify_one();
    }

    void pop(value_type* pItem) {
        boost::mutex::scoped_lock lock(mutex_);
        not_empty_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_empty, this));
        *pItem = container_[--unread_];
        lock.unlock();
        not_full_.notify_one();
    }

private:
    BoundedBuffer(const BoundedBuffer&);              // Disabled copy constructor
    BoundedBuffer& operator = (const BoundedBuffer&); // Disabled assign operator

    bool is_not_empty() const { return unread_ > 0; }
    bool is_not_full() const { return unread_ < container_.capacity(); }

    size_type unread_;
    container_type container_;
    boost::mutex mutex_;
    boost::condition not_empty_;
    boost::condition not_full_;
};
