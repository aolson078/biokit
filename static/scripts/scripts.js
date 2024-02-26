function showSubRadios(reportId) {
    var subRadios = document.querySelectorAll('.sub-radios');
    subRadios.forEach(function (subRadio) {
        subRadio.style.display = 'none';
    });
    var subRadioId = 'sub-' + reportId;
    var subRadioElement = document.getElementById(subRadioId);
    if (subRadioElement) {
        subRadioElement.style.display = 'block';
    }
}