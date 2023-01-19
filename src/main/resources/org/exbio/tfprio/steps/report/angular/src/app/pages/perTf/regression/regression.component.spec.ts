import { ComponentFixture, TestBed } from '@angular/core/testing';

import { RegressionComponent } from './regression.component';

describe('RegressionComponent', () => {
  let component: RegressionComponent;
  let fixture: ComponentFixture<RegressionComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ RegressionComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(RegressionComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
