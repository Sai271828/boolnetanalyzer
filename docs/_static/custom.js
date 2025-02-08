document.addEventListener('DOMContentLoaded', function() {
    // Create a search input element
    var searchInput = document.createElement('input');
    searchInput.setAttribute('type', 'text');
    searchInput.setAttribute('id', 'search-input');
    searchInput.setAttribute('placeholder', 'Search...');
    searchInput.style.margin = '10px';
    searchInput.style.padding = '5px';
    searchInput.style.width = '200px';

    // Create a search button element
    var searchButton = document.createElement('button');
    searchButton.innerText = 'Search';
    searchButton.style.margin = '10px';
    searchButton.style.padding = '5px';

    // Add the search input and button to the navigation bar
    var navBar = document.querySelector('nav');
    if (navBar) {
        navBar.appendChild(searchInput);
        navBar.appendChild(searchButton);
    }

    // Handle the search button click event
    searchButton.addEventListener('click', function() {
        var query = searchInput.value;
        if (query) {
            window.location.href = '/search.html?q=' + encodeURIComponent(query);
        }
    });

    // Handle enter key press in search input
    searchInput.addEventListener('keypress', function(e) {
        if (e.key === 'Enter') {
            var query = searchInput.value;
            if (query) {
                window.location.href = '/search.html?q=' + encodeURIComponent(query);
            }
        }
    });
});
