use std::sync::Mutex;

pub struct Pool<T: Default> {
    pool: Mutex<Vec<T>>,
}

impl<T: Default> Pool<T> {
    pub fn new() -> Self {
        Pool {
            pool: Mutex::new(Vec::new()),
        }
    }

    pub fn acquire(&self) -> PooledHandle<'_, T> {
        let value = self.pool.lock().unwrap().pop().unwrap_or_default();
        PooledHandle { value: Some(value), pool: self  }
    }

    pub fn release(&self, value: T) {
        self.pool.lock().unwrap().push(value);
    }

    /// Drop all currently retained items, freeing their memory.
    /// Items that are currently checked out are unaffected.
    pub fn clear(&self) {
        self.pool.lock().unwrap().clear();
    }
}

pub struct PooledHandle<'a, T: Default> {
    value: Option<T>,
    pool: &'a Pool<T>,
}

impl<'a, T: Default> std::ops::Deref for PooledHandle<'a, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        self.value.as_ref().unwrap()
    }
}

impl<'a, T: Default> std::ops::DerefMut for PooledHandle<'a, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.value.as_mut().unwrap()
    }
}

impl<'a, T: Default> Drop for PooledHandle<'a, T> {
    fn drop(&mut self) {
        if let Some(value) = self.value.take() {
            self.pool.release(value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn acquire_from_empty_pool_returns_default() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let handle = pool.acquire();
        assert!(handle.is_empty());
    }

    #[test]
    fn dropped_handle_returns_value_to_pool() {
        let pool: Pool<Vec<i32>> = Pool::new();
        {
            let mut handle = pool.acquire();
            handle.push(42);
        } // dropped here — should return to pool
        // Next acquire should get the same vec back (with its data)
        let handle = pool.acquire();
        assert_eq!(*handle, vec![42]);
    }

    #[test]
    fn deref_and_deref_mut_work() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let mut handle = pool.acquire();
        handle.push(1);
        handle.push(2);
        assert_eq!(handle.len(), 2);
        assert_eq!(handle[0], 1);
    }

    #[test]
    fn pool_reuses_allocations_across_sequential_acquires() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let ptr = {
            let mut h = pool.acquire();
            h.push(99);
            h.as_ptr() // pointer to the vec's heap allocation
        };
        // The same backing allocation should come back
        let h2 = pool.acquire();
        assert_eq!(h2.as_ptr(), ptr);
    }

    #[test]
    fn multiple_concurrent_acquires_get_distinct_handles() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let h1 = pool.acquire();
        let h2 = pool.acquire(); // pool is empty, gets a fresh default
        // Both handles are live simultaneously — they must be independent
        assert!(!std::ptr::eq(&*h1 as *const _, &*h2 as *const _));
    }

    #[test]
    fn pool_holds_multiple_returned_buffers() {
        let pool: Pool<Vec<i32>> = Pool::new();
        // Acquire two handles simultaneously, then drop both
        let h1 = pool.acquire();
        let h2 = pool.acquire();
        drop(h1);
        drop(h2);
        // Pool should now have two entries
        assert_eq!(pool.pool.lock().unwrap().len(), 2);
    }

    #[test]
    fn multithreaded_acquire_and_release() {
        let pool = Arc::new(Pool::<Vec<u64>>::new());
        let threads: Vec<_> = (0..8)
            .map(|i| {
                let pool = Arc::clone(&pool);
                std::thread::spawn(move || {
                    for _ in 0..100 {
                        let mut h = pool.acquire();
                        h.push(i as u64);
                        // handle dropped at end of loop iteration
                    }
                })
            })
            .collect();
        for t in threads {
            t.join().unwrap();
        }
        // All buffers should have been returned
        let n = pool.pool.lock().unwrap().len();
        assert!(n > 0, "expected buffers to be returned to pool, found {n}");
    }

    #[test]
    fn clear_empty_pool_is_a_no_op() {
        let pool: Pool<Vec<i32>> = Pool::new();
        pool.clear();
        assert_eq!(pool.pool.lock().unwrap().len(), 0);
    }

    #[test]
    fn clear_removes_retained_items() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let h1 = pool.acquire();
        let h2 = pool.acquire();
        drop(h1);
        drop(h2);
        assert_eq!(pool.pool.lock().unwrap().len(), 2);
        pool.clear();
        assert_eq!(pool.pool.lock().unwrap().len(), 0);
    }

    #[test]
    fn clear_does_not_affect_checked_out_handles() {
        let pool: Pool<Vec<i32>> = Pool::new();
        let mut h = pool.acquire();
        h.push(1);
        pool.clear(); // clears the pool but h is still live
        drop(h);      // h is returned after clear — pool now has 1 item
        assert_eq!(pool.pool.lock().unwrap().len(), 1);
    }

    #[test]
    fn acquire_after_clear_returns_default() {
        let pool: Pool<Vec<i32>> = Pool::new();
        {
            let mut h = pool.acquire();
            h.push(99);
        }
        pool.clear();
        let h = pool.acquire();
        assert!(h.is_empty(), "expected default after clear, got non-empty vec");
    }
}