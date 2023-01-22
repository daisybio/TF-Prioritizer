import { ComponentFixture, TestBed } from '@angular/core/testing';

import { RegressionCoefficientsComponent } from './regression-coefficients.component';

describe('RegressionCoefficientsComponent', () => {
  let component: RegressionCoefficientsComponent;
  let fixture: ComponentFixture<RegressionCoefficientsComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ RegressionCoefficientsComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(RegressionCoefficientsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
